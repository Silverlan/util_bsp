/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

module;

#include <fsys/filesystem.h>
#include <fsys/ifile.hpp>
#include <sharedutils/util_string.h>
#include <sharedutils/util_file.h>
#include <sharedutils/util_path.hpp>
#include <mathutil/uvec.h>
#ifdef ENABLE_VMT_SUPPORT
#include <VMTFile.h>
#endif
#include <vmf_entity_data.hpp>
#include <cassert>

module source_engine.bsp;

import util_zip;

#define IDBSPHEADER (('P' << 24) + ('S' << 16) + ('B' << 8) + 'V')

std::unique_ptr<source_engine::bsp::File> source_engine::bsp::File::Open(VFilePtr &f, ResultCode &code)
{
	if(f == nullptr) {
		code = ResultCode::FileNotFound;
		return nullptr;
	}
	auto header = f->Read<dheader_t>();
	if(header.ident != IDBSPHEADER) {
		code = ResultCode::InvalidHeaderIdent;
		return nullptr;
	}
	//if(header.version <= 19)
	//	throw std::runtime_error("BSP Map version " +std::to_string(header.version) +" is unsupported!");
	auto r = std::unique_ptr<File>(new File(f, header));
	code = ResultCode::Success;
	return r;
}

source_engine::bsp::File::File(VFilePtr &f, const dheader_t &header) : m_file(f), m_header(header) {}

bool source_engine::bsp::File::HasReadLump(uint32_t lumpId) const { return m_readLumps & (1ull << lumpId); }
void source_engine::bsp::File::MarkLumpAsRead(uint32_t lumpId) { m_readLumps |= (1ull << lumpId); }
#define LZMA_ID (('A' << 24) | ('M' << 16) | ('Z' << 8) | ('L'))
#pragma pack(push, 1)
struct lzma_header_t {
	uint32_t id;
	uint32_t actualSize; // always little endian
	uint32_t lzmaSize;   // always little endian
	std::array<uint8_t, 5> properties;
};
#pragma pack(pop)
enum class DecompressionResult : uint8_t { NotCompressed = 0, Success, Failed };
static DecompressionResult decompress_lzma(VFilePtr &f, std::vector<uint8_t> &decompressedData)
{
	auto id = f->Read<uint32_t>();
	f->Seek(f->Tell() - sizeof(uint32_t));
	if(id != LZMA_ID)
		return DecompressionResult::NotCompressed;
	auto lzmaHeader = f->Read<lzma_header_t>();
	std::vector<uint8_t> compressedData;
	compressedData.resize(lzmaHeader.lzmaSize);
	f->Read(compressedData.data(), compressedData.size());

	decompressedData.resize(lzmaHeader.actualSize);
	size_t decompressedSize = decompressedData.size();
	size_t compressedSize = compressedData.size();
	auto result = source_engine::bsp::lzma_uncompress(decompressedData.data(), &decompressedSize, compressedData.data(), &compressedSize, lzmaHeader.properties.data(), lzmaHeader.properties.size());
	return result ? DecompressionResult::Success : DecompressionResult::Failed;
}

template<class TLump>
static std::unique_ptr<ufile::IFile> get_lump_file(const TLump &lump, VFilePtr &f, uint64_t &offset, uint64_t &outSize)
{
	std::vector<uint8_t> uncompressedData;
	auto r = decompress_lzma(f, uncompressedData);
	if(r == DecompressionResult::Success) {
		offset = 0;
		outSize = uncompressedData.size();
		return std::make_unique<ufile::VectorFile>(std::move(uncompressedData));
	}
	if(r == DecompressionResult::Failed)
		return nullptr;
	offset = lump.fileofs;
	outSize = lump.filelen;
	return std::make_unique<fsys::File>(f);
}

template<class T, class TContainer>
void source_engine::bsp::File::ReadData(uint32_t lumpId, TContainer &data)
{
	if(HasReadLump(lumpId) == true)
		return;
	MarkLumpAsRead(lumpId);
	auto &f = m_file;
	auto &header = m_header;
	auto &lump = header.lumps.at(lumpId);
	if(lump.filelen == 0)
		return;
	f->Seek(lump.fileofs);

	size_t offset;
	size_t size;
	auto fl = get_lump_file(lump, f, offset, size);
	if(!fl)
		return;

	assert((size % sizeof(T)) == 0);
	auto numData = size / sizeof(T);
	data.resize(numData);
	fl->Read(&data[0], sizeof(data.front()) * data.size());
}

template<class T, class TContainer>
void source_engine::bsp::File::ReadData(uint32_t lumpId, TContainer &data, uint32_t padding)
{
	if(padding == 0) {
		ReadData<T, TContainer>(lumpId, data);
		return;
	}
	if(HasReadLump(lumpId) == true)
		return;
	MarkLumpAsRead(lumpId);
	auto &f = m_file;
	auto &header = m_header;
	auto &lump = header.lumps.at(lumpId);
	if(lump.filelen == 0)
		return;
	f->Seek(lump.fileofs);

	size_t offset;
	size_t size;
	auto fl = get_lump_file(lump, f, offset, size);
	if(!fl)
		return;

	auto id = fl->Read<uint32_t>();
	fl->Seek(fl->Tell() - sizeof(uint32_t));

	assert((size % sizeof(T)) == 0);
	auto numData = size / sizeof(T);
	data.resize(numData);
	for(auto i = decltype(numData) {0}; i < numData; ++i) {
		fl->Read(&data[i], sizeof(data.front()));
		fl->Seek(fl->Tell() + padding);
	}
}

void source_engine::bsp::File::ReadGameData()
{
	const auto lumpId = LUMP_ID_GAME;
	if(HasReadLump(lumpId) == true)
		return;
	MarkLumpAsRead(lumpId);
	auto &f = m_file;
	auto &header = m_header;
	auto &lump = header.lumps.at(lumpId);
	if(lump.filelen == 0)
		return;
	f->Seek(lump.fileofs);

	auto numLumps = f->Read<int32_t>();
	m_gameLumps.resize(numLumps);

	for(auto i = decltype(numLumps) {0u}; i < numLumps; ++i) {
		auto &lump = m_gameLumps[i];
		size_t offset;
		size_t size;
		auto fl = get_lump_file(lump, f, offset, size);
		if(!fl)
			return;
		fl->Read(&lump, sizeof(lump));
	}
}

void source_engine::bsp::File::ReadStaticPropsData()
{
	ReadGameData();
	auto it = std::find_if(m_gameLumps.begin(), m_gameLumps.end(), [](const dgamelump_t &lump) {
		std::string identifier(reinterpret_cast<const char *>(&lump.id), 4);
		std::reverse(identifier.begin(), identifier.end());
		return identifier == "sprp";
	});
	if(it == m_gameLumps.end())
		return;
	auto &f = m_file;
	auto &lump = *it;
	if(lump.filelen == 0)
		return;
	auto version = lump.version;
	f->Seek(lump.fileofs);

	size_t offset;
	size_t size;
	auto fl = get_lump_file(lump, f, offset, size);
	if(!fl)
		return;
	auto numDictEntires = fl->Read<int32_t>();
	auto &modelNames = m_staticPropData.dictionaryModelNames;
	modelNames.reserve(numDictEntires);
	for(auto i = decltype(numDictEntires) {0}; i < numDictEntires; ++i) {
		std::array<char, 128> name;
		fl->Read(name.data(), name.size() * sizeof(name.front()));
		modelNames.push_back(std::string {name.data()});
	}

	auto numLeafEntries = fl->Read<int32_t>();
	auto &leaves = m_staticPropData.leaves;
	leaves.resize(numLeafEntries);
	fl->Read(leaves.data(), leaves.size() * sizeof(leaves.front()));

	auto numStaticProps = fl->Read<int32_t>();
	auto szStaticProps = size - (fl->Tell() - offset);
	assert((szStaticProps % numStaticProps) == 0);
	auto szPerLump = (numStaticProps > 0) ? (szStaticProps / numStaticProps) : 0;
	auto &staticPropLumps = m_staticPropData.staticPropLumps;
	staticPropLumps.reserve(numStaticProps);
	for(auto i = decltype(numStaticProps) {0u}; i < numStaticProps; ++i) {
		staticPropLumps.push_back({});
		auto offset = fl->Tell();
		auto &propLump = staticPropLumps.back();
		if(version >= 4) {
			propLump.Origin = fl->Read<Vector3>();
			propLump.Angles = fl->Read<EulerAngles>();
			propLump.PropType = fl->Read<uint16_t>();
			propLump.FirstLeaf = fl->Read<uint16_t>();
			propLump.LeafCount = fl->Read<uint16_t>();
			propLump.Solid = fl->Read<uint8_t>();
			propLump.Flags = fl->Read<uint8_t>();
			propLump.Skin = fl->Read<int32_t>();
			propLump.FadeMinDist = fl->Read<float>();
			propLump.FadeMaxDist = fl->Read<float>();
			propLump.LightingOrigin = fl->Read<Vector3>();
		}
		fl->Seek(offset + szPerLump);
	}
}

void source_engine::bsp::File::ReadEntityData()
{
	const auto lumpId = LUMP_ID_ENTITIES;
	if(HasReadLump(lumpId) == true)
		return;
	MarkLumpAsRead(lumpId);
	auto &f = m_file;
	auto &header = m_header;
	auto &lump = header.lumps.at(lumpId);
	if(lump.filelen == 0)
		return;
	f->Seek(lump.fileofs);

	size_t offset;
	size_t size;
	auto fl = get_lump_file(lump, f, offset, size);
	if(!fl)
		return;

	auto id = fl->Read<uint32_t>();
	fl->Seek(fl->Tell() - sizeof(uint32_t));

	auto lumpEnd = offset + size;
	std::unique_ptr<vmf::DataFileBlock> data {vmf::DataFile::ReadBlock(*fl, lumpEnd)};
	if(data == nullptr)
		return;
	auto it = data->blocks.find("unnamed");
	if(it == data->blocks.end() || it->second == nullptr)
		return;
	m_entities.reserve(it->second->size());
	for(auto *block : *it->second)
		m_entities.push_back(EntityBlock(block));
	it->second->clear();

	//m_entities.resize(lump.filelen);
	//f->Read(&m_entities[0],lump.filelen);
}

void source_engine::bsp::File::ReadLightMapData()
{
	const auto lumpId = LUMP_ID_LIGHTING;
	if(HasReadLump(lumpId) == true)
		return;
	MarkLumpAsRead(lumpId);

	auto &f = m_file;
	auto &header = m_header;
	auto &lump = header.lumps.at(lumpId);
	if(lump.filelen == 0)
		return;
	f->Seek(lump.fileofs);

	size_t offset;
	size_t size;
	auto fl = get_lump_file(lump, f, offset, size);
	if(!fl)
		return;
	m_lightMapData.resize(size);
	fl->Read(m_lightMapData.data(), size);
}
void source_engine::bsp::File::ReadHDRLightMapData()
{
	const auto lumpId = LUMP_ID_LIGHTING_HDR;
	if(HasReadLump(lumpId) == true)
		return;
	MarkLumpAsRead(lumpId);

	auto &f = m_file;
	auto &header = m_header;
	auto &lump = header.lumps.at(lumpId);
	if(lump.filelen == 0)
		return;
	f->Seek(lump.fileofs);

	size_t offset;
	size_t size;
	auto fl = get_lump_file(lump, f, offset, size);
	if(!fl)
		return;
	m_lightMapDataHDR.resize(size);
	fl->Read(m_lightMapDataHDR.data(), size);
}
void source_engine::bsp::File::ReadVisibilityData()
{
	const auto lumpId = LUMP_ID_VISIBILITY;
	if(HasReadLump(lumpId) == true)
		return;
	MarkLumpAsRead(lumpId);

	auto &f = m_file;
	auto &header = m_header;
	auto &lump = header.lumps.at(lumpId);
	if(lump.filelen == 0)
		return;
	f->Seek(lump.fileofs);

	size_t offset;
	size_t size;
	auto fl = get_lump_file(lump, f, offset, size);
	if(!fl)
		return;

#pragma pack(push, 1)
	struct ClusterOffsetInfo {
		int32_t offsetPVS;
		int32_t offsetPAS;
	};
#pragma pack(pop)
	auto numClusters = fl->Read<int32_t>();
	std::vector<ClusterOffsetInfo> offsets(numClusters);
	fl->Read(offsets.data(), offsets.size() * sizeof(offsets.front()));

	auto headerSize = fl->Tell() - offset;
	std::vector<uint8_t> compressedData(size - headerSize);
	fl->Read(compressedData.data(), compressedData.size() * sizeof(compressedData.front()));
	m_visibilityData.resize(numClusters, std::vector<uint8_t>(numClusters, 0u));
	for(auto i = decltype(numClusters) {0}; i < numClusters; ++i) {
		auto v = offsets.at(i).offsetPVS - headerSize;
		auto &clusterVisible = m_visibilityData.at(i);
		for(auto c = decltype(numClusters) {0u}; c < numClusters; v++) {
			if(compressedData.at(v) == 0) {
				v++;
				c += 8 * compressedData.at(v);
			}
			else {
				for(uint8_t bit = 1; bit != 0; bit *= 2, c++) {
					if(compressedData.at(v) & bit)
						clusterVisible.at(c) = 1;
				}
			}
		}
	}
}
void source_engine::bsp::File::ReadPlaneData() { ReadData<std::remove_reference_t<decltype(m_planes.front())>>(LUMP_ID_PLANES, m_planes); }
void source_engine::bsp::File::ReadVertexData() { ReadData<std::remove_reference_t<decltype(m_vertices.front())>>(LUMP_ID_VERTICES, m_vertices); }
void source_engine::bsp::File::ReadEdgeData() { ReadData<std::remove_reference_t<decltype(m_edges.front())>>(LUMP_ID_EDGES, m_edges); }
void source_engine::bsp::File::ReadSurfEdgeData() { ReadData<std::remove_reference_t<decltype(m_surfEdges.front())>>(LUMP_ID_SURF_EDGES, m_surfEdges); }
void source_engine::bsp::File::ReadFaceData() { ReadData<std::remove_reference_t<decltype(m_faces.front())>>(LUMP_ID_FACES, m_faces); }
void source_engine::bsp::File::ReadHDRFaceData() { ReadData<std::remove_reference_t<decltype(m_hdrFaces.front())>>(LUMP_ID_FACES_HDR, m_hdrFaces); }
void source_engine::bsp::File::ReadOriginalFaceData() { ReadData<std::remove_reference_t<decltype(m_origFaces.front())>>(LUMP_ID_ORIGINAL_FACES, m_origFaces); }
void source_engine::bsp::File::ReadBrushData() { ReadData<std::remove_reference_t<decltype(m_brushes.front())>>(LUMP_ID_BRUSHES, m_brushes); }
void source_engine::bsp::File::ReadBrushSideData() { ReadData<std::remove_reference_t<decltype(m_brushSides.front())>>(LUMP_ID_BRUSH_SIDES, m_brushSides); }
void source_engine::bsp::File::ReadTexInfoData() { ReadData<std::remove_reference_t<decltype(m_texInfo.front())>>(LUMP_ID_TEXINFO, m_texInfo); }
void source_engine::bsp::File::ReadTexData() { ReadData<std::remove_reference_t<decltype(m_texData.front())>>(LUMP_ID_TEXDATA, m_texData); }
void source_engine::bsp::File::ReadModelData() { ReadData<std::remove_reference_t<decltype(m_models.front())>>(LUMP_ID_MODELS, m_models); }
void source_engine::bsp::File::ReadDispLightmapSamplePositions() { ReadData<std::remove_reference_t<decltype(m_dispLightmapSamplePositions.front())>>(LUMP_ID_DISP_LIGHTMAP_SAMPLE_POSITIONS, m_dispLightmapSamplePositions); }
void source_engine::bsp::File::ReadCubemapSamples() { ReadData<std::remove_reference_t<decltype(m_cubemapSamples.front())>>(LUMP_ID_CUBEMAPS, m_cubemapSamples); }
void source_engine::bsp::File::ReadTexDataStringTableData() { ReadData<std::remove_reference_t<decltype(m_texDataStringTableData.front())>>(LUMP_ID_TEXDATA_STRING_TABLE, m_texDataStringTableData); }
void source_engine::bsp::File::ReadTexDataStringData()
{
	const auto lumpId = LUMP_ID_TEXDATA_STRING_DATA;
	if(HasReadLump(lumpId) == true)
		return;
	MarkLumpAsRead(lumpId);
	auto &f = m_file;
	auto &header = m_header;
	auto &lump = header.lumps.at(lumpId);
	if(lump.filelen == 0)
		return;
	f->Seek(lump.fileofs);

	size_t offset;
	size_t size;
	auto fl = get_lump_file(lump, f, offset, size);
	if(!fl)
		return;

	while(fl->Eof() == false && fl->Tell() < offset + size) {
		m_texDataStringDataIndexMap.insert(std::make_pair(fl->Tell() - offset, m_texDataStringData.size()));
		m_texDataStringData.push_back(::util::Path::CreateFile(fl->ReadString()).GetString());
	}
}
void source_engine::bsp::File::ReadNodes() { ReadData<std::remove_reference_t<decltype(m_nodes.front())>>(5u, m_nodes); }
void source_engine::bsp::File::ReadLeaves()
{
	auto padding = 0u;
	if(m_header.version <= 19)
		padding = 56u - sizeof(decltype(m_leaves.front()));
	ReadData<std::remove_reference_t<decltype(m_leaves.front())>>(10u, m_leaves, padding);
}
void source_engine::bsp::File::ReadLeafFaces() { ReadData<std::remove_reference_t<decltype(m_leafFaces.front())>>(16u, m_leafFaces); }
void source_engine::bsp::File::ReadLeafBrushes() { ReadData<std::remove_reference_t<decltype(m_leafBrushes.front())>>(17u, m_leafBrushes); }
#define PKID(a, b) (((b) << 24) | ((a) << 16) | ('K' << 8) | 'P')

void source_engine::bsp::File::ReadPakfile()
{
	const auto lumpId = LUMP_ID_PAKFILE;
	if(HasReadLump(lumpId) == true)
		return;
	MarkLumpAsRead(lumpId);
	auto &f = m_file;
	auto &header = m_header;
	auto &lump = header.lumps.at(lumpId);
	if(lump.filelen == 0)
		return;

	int64_t offset = lump.fileofs + lump.filelen - sizeof(ZIP_EndOfCentralDirRecord);

	auto &data = m_pakZipData;
	m_file->Seek(lump.fileofs);
	data.resize(lump.filelen);
	f->Read(data.data(), data.size());
	m_pakZipFile = uzip::ZIPFile::Open(data.data(), data.size());

	auto bFoundRecord = false;
	for(; offset >= 0; offset--) {
		f->Seek(offset);
		f->Read(&m_zipDirRecord, sizeof(m_zipDirRecord));
		if(m_zipDirRecord.signature == PKID(5, 6)) {
			bFoundRecord = true;
			break;
		}
	}
	if(bFoundRecord == false)
		return;
	f->Seek(lump.fileofs + m_zipDirRecord.startOfCentralDirOffset);
	m_fileHeaders.resize(m_zipDirRecord.nCentralDirectoryEntries_Total);
	m_localFileHeaders.resize(m_fileHeaders.size());
	m_fileNames.reserve(m_fileHeaders.size());
	m_fileDataOffsets.reserve(m_fileHeaders.size());
	for(auto i = decltype(m_fileHeaders.size()) {0}; i < m_fileHeaders.size(); ++i) {
		auto &hd = m_fileHeaders.at(i);
		auto &lhd = m_localFileHeaders.at(i);
		auto offset = f->Tell();
		f->Read(&hd, sizeof(hd));

		m_fileNames.push_back({});
		auto &fname = m_fileNames.back();
		fname.resize(hd.fileNameLength);
		f->Read(&fname[0], hd.fileNameLength);

		f->Seek(lump.fileofs + hd.relativeOffsetOfLocalHeader);
		f->Read(&lhd, sizeof(lhd));
		m_fileDataOffsets.push_back(f->Tell() + hd.fileNameLength + hd.extraFieldLength);

		f->Seek(offset + sizeof(hd) + hd.fileNameLength + hd.extraFieldLength + hd.fileCommentLength);
	}
}
void source_engine::bsp::File::ReadDisplacementData()
{
	if(m_bDisplacementsRead == true)
		return;
	m_bDisplacementsRead = true;
	auto &displacements = m_displacements;
	auto &faces = GetFaces();
	auto &origFaces = GetOriginalFaces();
	auto &dispInfo = GetDispInfo();
	auto &planes = GetPlanes();
	for(auto &face : origFaces) {
		if(face.dispinfo == -1)
			continue;
		auto it = std::find_if(displacements.begin(), displacements.end(), [&face](const dDisp &disp) { return (disp.face.dispinfo == face.dispinfo) ? true : false; });
		if(it != displacements.end())
			continue;
		auto &plane = planes.at(face.planenum);
		auto &faceDispInfo = dispInfo.at(face.dispinfo);
		if(displacements.size() == displacements.capacity())
			displacements.reserve(displacements.size() + 100);
		displacements.push_back({faces.at(faceDispInfo.MapFace), plane, faceDispInfo});
	}

	auto &header = m_header;
	size_t offsetV;
	size_t sizeV;
	auto &lumpV = header.lumps.at(LUMP_ID_DISP_VERTS);
	m_file->Seek(lumpV.fileofs);
	auto flV = get_lump_file(lumpV, m_file, offsetV, sizeV);

	size_t offsetT;
	size_t sizeT;
	auto &lumpT = header.lumps.at(LUMP_ID_DISP_TRIS);
	m_file->Seek(lumpT.fileofs);
	auto flT = get_lump_file(lumpT, m_file, offsetT, sizeT);
	if(!flV || !flT)
		return;
	for(auto &disp : displacements) {
		auto &dispInfo = disp.dispInfo;
		auto &dispVerts = disp.verts;
		auto &dispTris = disp.tris;
		auto p = static_cast<uint32_t>(dispInfo.power);
		auto numVerts = umath::pow2(umath::pow(2u, p) + 1u);
		auto numTris = 2 * umath::pow2(umath::pow(2u, p));

		if(sizeV != 0) {
			dispVerts.resize(numVerts);
			flV->Seek(offsetV + dispInfo.DispVertStart * sizeof(dispVerts.front()));
			flV->Read(dispVerts.data(), dispVerts.size() * sizeof(dispVerts.front()));
		}

		if(sizeT != 0) {
			dispTris.resize(numTris);
			flT->Seek(offsetT + dispInfo.DispTriStart * sizeof(dispTris.front()));
			flT->Read(dispTris.data(), dispTris.size() * sizeof(dispTris.front()));
		}
	}
}
void source_engine::bsp::File::ReadDispInfo() { ReadData<std::remove_reference_t<decltype(m_dispInfo.front())>>(26u, m_dispInfo); }

const source_engine::bsp::lump_t *source_engine::bsp::File::GetLumpHeaderInfo(uint32_t lumpId) const
{
	auto &header = m_header;
	return (lumpId < header.lumps.size()) ? &header.lumps.at(lumpId) : nullptr;
}
const std::vector<source_engine::bsp::dgamelump_t> &source_engine::bsp::File::GetGameLumps()
{
	ReadGameData();
	return m_gameLumps;
}
const std::vector<source_engine::bsp::EntityBlock> &source_engine::bsp::File::GetEntities()
{
	ReadEntityData();
	return m_entities;
}
const std::vector<source_engine::bsp::dplane_t> &source_engine::bsp::File::GetPlanes()
{
	ReadPlaneData();
	return m_planes;
}
const std::vector<Vector3> &source_engine::bsp::File::GetVertices()
{
	ReadVertexData();
	return m_vertices;
}
const std::vector<source_engine::bsp::dedge_t> &source_engine::bsp::File::GetEdges()
{
	ReadEdgeData();
	return m_edges;
}
const std::vector<int32_t> &source_engine::bsp::File::GetSurfEdges()
{
	ReadSurfEdgeData();
	return m_surfEdges;
}
const std::vector<source_engine::bsp::dface_t> &source_engine::bsp::File::GetFaces()
{
	ReadFaceData();
	return m_faces;
}
const std::vector<source_engine::bsp::dface_t> &source_engine::bsp::File::GetHDRFaces()
{
	ReadHDRFaceData();
	return m_hdrFaces;
}
const std::vector<source_engine::bsp::dface_t> &source_engine::bsp::File::GetOriginalFaces()
{
	ReadOriginalFaceData();
	return m_origFaces;
}
const std::vector<source_engine::bsp::dbrush_t> &source_engine::bsp::File::GetBrushes()
{
	ReadBrushData();
	return m_brushes;
}
const std::vector<source_engine::bsp::dbrushside_t> &source_engine::bsp::File::GetBrushSides()
{
	ReadBrushSideData();
	return m_brushSides;
}
const std::vector<source_engine::bsp::texinfo_t> &source_engine::bsp::File::GetTexInfo()
{
	ReadTexInfoData();
	return m_texInfo;
}
const std::vector<source_engine::bsp::dtexdata_t> &source_engine::bsp::File::GetTexData()
{
	ReadTexData();
	return m_texData;
}
const std::vector<source_engine::bsp::dmodel_t> &source_engine::bsp::File::GetModels()
{
	ReadModelData();
	return m_models;
}
const std::vector<std::string> &source_engine::bsp::File::GetTexDataStrings()
{
	ReadTexDataStringData();
	return m_texDataStringData;
}
const std::vector<std::string> &source_engine::bsp::File::GetTranslatedTexDataStrings()
{
	if(m_bStringDataTranslated == true)
		return m_texDataStringDataTranslated;
	m_bStringDataTranslated = true;
	auto &strings = m_texDataStringDataTranslated = GetTexDataStrings();
#if ENABLE_VMT_SUPPORT
	auto &fileNames = GetFilenames();
	for(auto &fname : fileNames) {
		auto relFname = fname;
		ufile::remove_extension_from_filename(relFname);
		if(ustring::compare(relFname.c_str(), "materials/", false, 10))
			relFname = ustring::substr(relFname, 10);

		auto it = std::find_if(strings.begin(), strings.end(), [&relFname](const std::string &fnameOther) { return ustring::compare(relFname, fnameOther, false); });
		if(it == strings.end())
			continue;
		std::string ext;
		if(ufile::get_extension(fname, &ext) == true && ustring::compare<std::string>(ext, "vmt", false) == true) {
			std::vector<uint8_t> data;
			if(ReadFile(fname, data) == true) {
				VTFLib::CVMTFile vmt {};
				if(vmt.Load(data.data(), static_cast<vlUInt>(data.size())) == vlTrue) {
					auto *vmtRoot = vmt.GetRoot();
					std::string shader = vmtRoot->GetName();
					if(ustring::compare<std::string>(shader, "patch", false)) {
						VTFLib::Nodes::CVMTNode *node = nullptr;
						if((node = vmtRoot->GetNode("include")) != nullptr) {
							if(node->GetType() == VMTNodeType::NODE_TYPE_STRING) {
								auto *includeNode = static_cast<VTFLib::Nodes::CVMTStringNode *>(node);
								std::string include = includeNode->GetValue();
								ufile::remove_extension_from_filename(include);
								if(ustring::compare(include.c_str(), "materials/", false, 10))
									*it = ustring::substr(include, 10);
							}
						}
					}
				}
			}
		}
	}
#endif
	return strings;
}
const std::vector<uint32_t> &source_engine::bsp::File::GetTexDataStringIndices()
{
	if(m_bStringDataIndicesRead == true)
		return m_texDataStringDataIndices;
	m_bStringDataIndicesRead = true;

	ReadTexDataStringTableData();
	ReadTexDataStringData();

	m_texDataStringDataIndices.reserve(m_texDataStringDataIndexMap.size());
	for(auto offset : m_texDataStringTableData) {
		auto it = m_texDataStringDataIndexMap.find(offset);
		assert(it != m_texDataStringDataIndexMap.end());
		if(it == m_texDataStringDataIndexMap.end())
			m_texDataStringDataIndices.push_back(0u);
		else
			m_texDataStringDataIndices.push_back(it->second);
	}
	return m_texDataStringDataIndices;
}
const std::vector<source_engine::bsp::dnode_t> &source_engine::bsp::File::GetNodes()
{
	ReadNodes();
	return m_nodes;
}
const std::vector<source_engine::bsp::dleaf_t> &source_engine::bsp::File::GetLeaves()
{
	ReadLeaves();
	return m_leaves;
}
const std::vector<uint16_t> &source_engine::bsp::File::GetLeafFaces()
{
	ReadLeafFaces();
	return m_leafFaces;
}
const std::vector<uint16_t> &source_engine::bsp::File::GetLeafBrushes()
{
	ReadLeafBrushes();
	return m_leafBrushes;
}
const std::vector<std::string> &source_engine::bsp::File::GetFilenames()
{
	ReadPakfile();
	return m_fileNames;
}
const std::vector<source_engine::bsp::ddispinfo_t> &source_engine::bsp::File::GetDispInfo()
{
	ReadDispInfo();
	return m_dispInfo;
}
const std::vector<source_engine::bsp::dDisp> &source_engine::bsp::File::GetDisplacements()
{
	ReadDisplacementData();
	return m_displacements;
}
const std::vector<uint8_t> &source_engine::bsp::File::GetLightMapData()
{
	ReadLightMapData();
	return m_lightMapData;
}
const std::vector<uint8_t> &source_engine::bsp::File::GetHDRLightMapData()
{
	ReadHDRLightMapData();
	return m_lightMapDataHDR;
}
const std::vector<std::vector<uint8_t>> &source_engine::bsp::File::GetVisibilityData()
{
	ReadVisibilityData();
	return m_visibilityData;
}
const std::vector<uint8_t> &source_engine::bsp::File::GetDispLightmapSamplePositions()
{
	ReadDispLightmapSamplePositions();
	return m_dispLightmapSamplePositions;
}
const std::vector<source_engine::bsp::dcubemapsample_t> &source_engine::bsp::File::GetCubemapSamples()
{
	ReadCubemapSamples();
	return m_cubemapSamples;
}
const source_engine::bsp::StaticPropData &source_engine::bsp::File::GetStaticPropData()
{
	ReadStaticPropsData();
	return m_staticPropData;
}
const source_engine::bsp::dheader_t &source_engine::bsp::File::GetHeaderData() const { return m_header; }
const bool source_engine::bsp::File::ReadFile(const std::string &fname, std::vector<uint8_t> &data)
{
	ReadPakfile();
	if(!m_pakZipFile)
		return false;
	std::string err;
	auto r = m_pakZipFile->ReadFile(fname, data, err);
	if(r == false) {
		// Probably an unsupported compression algorithm
		// std::cout<<"Unable to read zip file '"<<fname<<"': "<<err<<std::endl;
	}
	return r;

	// Obsolete
#if 0
	auto it = std::find_if(m_fileNames.begin(),m_fileNames.end(),[&](const std::string &fnameOther) {
		return ustring::compare(fname,fnameOther,false);
	});
	if(it == m_fileNames.end())
		return false;
	auto hdIdx = it -m_fileNames.begin();
	auto &lhd = m_fileHeaders.at(hdIdx);
	auto &f = m_file;
	if(lhd.compressionMethod != 0 || lhd.uncompressedSize != lhd.compressedSize)
		return false;
	data.resize(lhd.uncompressedSize);
	f->Seek(m_fileDataOffsets.at(hdIdx));
	f->Read(data.data(),data.size());
	return true;
#endif
}
