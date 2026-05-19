// SPDX-FileCopyrightText: (c) 2024 Silverlan <opensource@pragma-engine.com>
// SPDX-License-Identifier: MIT

module;

#ifdef ENABLE_VMT_SUPPORT
#include <VMTFile.h>
#endif
#include <cassert>

module source_engine.bsp;

import util_zip;
import source_engine.vmf;

std::unique_ptr<source_engine::bsp::File> source_engine::bsp::File::Open(pragma::fs::VFilePtr &f, ResultCode &code)
{
	if(f == nullptr) {
		code = ResultCode::FileNotFound;
		return nullptr;
	}
	auto header = f->Read<Header>();
	constexpr int32_t expectedHeader = (('P' << 24) + ('S' << 16) + ('B' << 8) + 'V');
	if(header.identifier != expectedHeader) {
		code = ResultCode::InvalidHeaderIdent;
		return nullptr;
	}
	//if(header.version <= 19)
	//	throw std::runtime_error("BSP Map version " +std::to_string(header.version) +" is unsupported!");
	auto r = std::unique_ptr<File>(new File(f, header));
	code = ResultCode::Success;
	return r;
}

source_engine::bsp::File::File(pragma::fs::VFilePtr &f, const Header &header) : m_file(f), m_header(header) {}

bool source_engine::bsp::File::HasReadLump(LumpId lumpId) const { return m_readLumps & (1ull << pragma::math::to_integral(lumpId)); }
void source_engine::bsp::File::MarkLumpAsRead(LumpId lumpId) { m_readLumps |= (1ull << pragma::math::to_integral(lumpId)); }

#pragma pack(push, 1)
struct LzmaHeader {
	uint32_t id;
	uint32_t actualSize;
	uint32_t lzmaSize;
	std::array<uint8_t, 5> properties;
};
#pragma pack(pop)
enum class DecompressionResult : uint8_t { NotCompressed = 0, Success, Failed };
static DecompressionResult decompress_lzma(pragma::fs::VFilePtr &f, std::vector<uint8_t> &decompressedData)
{
	auto id = f->Read<uint32_t>();
	f->Seek(f->Tell() - sizeof(uint32_t));
	constexpr uint32_t expectedId = (('A' << 24) | ('M' << 16) | ('Z' << 8) | ('L'));
	if(id != expectedId)
		return DecompressionResult::NotCompressed;
	auto lzmaHeader = f->Read<LzmaHeader>();
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
static std::unique_ptr<ufile::IFile> get_lump_file(const TLump &lump, pragma::fs::VFilePtr &f, uint64_t &offset, uint64_t &outSize)
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
	offset = lump.fileOffset;
	outSize = lump.fileLength;
	return std::make_unique<pragma::fs::File>(f);
}

template<class T, class TContainer>
void source_engine::bsp::File::ReadData(LumpId lumpId, TContainer &data)
{
	if(HasReadLump(lumpId) == true)
		return;
	MarkLumpAsRead(lumpId);
	auto &f = m_file;
	auto &header = m_header;
	auto &lump = header.lumps.at(pragma::math::to_integral(lumpId));
	if(lump.fileLength == 0)
		return;
	f->Seek(lump.fileOffset);

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
void source_engine::bsp::File::ReadData(LumpId lumpId, TContainer &data, uint32_t padding)
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
	auto &lump = header.lumps.at(pragma::math::to_integral(lumpId));
	if(lump.fileLength == 0)
		return;
	f->Seek(lump.fileOffset);

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
	const auto lumpId = LumpId::Game;
	if(HasReadLump(lumpId) == true)
		return;
	MarkLumpAsRead(lumpId);
	auto &f = m_file;
	auto &header = m_header;
	auto &lump = header.lumps.at(pragma::math::to_integral(lumpId));
	if(lump.fileLength == 0)
		return;
	f->Seek(lump.fileOffset);

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
	auto it = std::find_if(m_gameLumps.begin(), m_gameLumps.end(), [](const GameLump &lump) {
		std::string identifier(reinterpret_cast<const char *>(&lump.id), 4);
		std::reverse(identifier.begin(), identifier.end());
		return identifier == "sprp";
	});
	if(it == m_gameLumps.end())
		return;
	auto &f = m_file;
	auto &lump = *it;
	if(lump.fileLength == 0)
		return;
	auto version = lump.version;
	f->Seek(lump.fileOffset);

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
			propLump.origin = fl->Read<Vector3>();
			propLump.angles = fl->Read<EulerAngles>();
			propLump.propType = fl->Read<uint16_t>();
			propLump.firstLeaf = fl->Read<uint16_t>();
			propLump.leafCount = fl->Read<uint16_t>();
			propLump.solid = fl->Read<uint8_t>();
			propLump.flags = fl->Read<uint8_t>();
			propLump.skin = fl->Read<int32_t>();
			propLump.fadeMinDist = fl->Read<float>();
			propLump.fadeMaxDist = fl->Read<float>();
			propLump.lightingOrigin = fl->Read<Vector3>();
		}
		fl->Seek(offset + szPerLump);
	}
}

void source_engine::bsp::File::ReadEntityData()
{
	const auto lumpId = LumpId::Entities;
	if(HasReadLump(lumpId) == true)
		return;
	MarkLumpAsRead(lumpId);
	auto &f = m_file;
	auto &header = m_header;
	auto &lump = header.lumps.at(pragma::math::to_integral(lumpId));
	if(lump.fileLength == 0)
		return;
	f->Seek(lump.fileOffset);

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

	//m_entities.resize(lump.fileLength);
	//f->Read(&m_entities[0],lump.fileLength);
}

void source_engine::bsp::File::ReadLightMapData()
{
	const auto lumpId = LumpId::Lighting;
	if(HasReadLump(lumpId) == true)
		return;
	MarkLumpAsRead(lumpId);

	auto &f = m_file;
	auto &header = m_header;
	auto &lump = header.lumps.at(pragma::math::to_integral(lumpId));
	if(lump.fileLength == 0)
		return;
	f->Seek(lump.fileOffset);

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
	const auto lumpId = LumpId::LightingHdr;
	if(HasReadLump(lumpId) == true)
		return;
	MarkLumpAsRead(lumpId);

	auto &f = m_file;
	auto &header = m_header;
	auto &lump = header.lumps.at(pragma::math::to_integral(lumpId));
	if(lump.fileLength == 0)
		return;
	f->Seek(lump.fileOffset);

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
	const auto lumpId = LumpId::Visibility;
	if(HasReadLump(lumpId) == true)
		return;
	MarkLumpAsRead(lumpId);

	auto &f = m_file;
	auto &header = m_header;
	auto &lump = header.lumps.at(pragma::math::to_integral(lumpId));
	if(lump.fileLength == 0)
		return;
	f->Seek(lump.fileOffset);

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
void source_engine::bsp::File::ReadPlaneData() { ReadData<std::remove_reference_t<decltype(m_planes.front())>>(LumpId::Planes, m_planes); }
void source_engine::bsp::File::ReadVertexData() { ReadData<std::remove_reference_t<decltype(m_vertices.front())>>(LumpId::Vertices, m_vertices); }
void source_engine::bsp::File::ReadEdgeData() { ReadData<std::remove_reference_t<decltype(m_edges.front())>>(LumpId::Edges, m_edges); }
void source_engine::bsp::File::ReadSurfEdgeData() { ReadData<std::remove_reference_t<decltype(m_surfEdges.front())>>(LumpId::SurfEdges, m_surfEdges); }
void source_engine::bsp::File::ReadFaceData() { ReadData<std::remove_reference_t<decltype(m_faces.front())>>(LumpId::Faces, m_faces); }
void source_engine::bsp::File::ReadHDRFaceData() { ReadData<std::remove_reference_t<decltype(m_hdrFaces.front())>>(LumpId::FacesHdr, m_hdrFaces); }
void source_engine::bsp::File::ReadOriginalFaceData() { ReadData<std::remove_reference_t<decltype(m_origFaces.front())>>(LumpId::OriginalFaces, m_origFaces); }
void source_engine::bsp::File::ReadBrushData() { ReadData<std::remove_reference_t<decltype(m_brushes.front())>>(LumpId::Brushes, m_brushes); }
void source_engine::bsp::File::ReadBrushSideData() { ReadData<std::remove_reference_t<decltype(m_brushSides.front())>>(LumpId::BrushSides, m_brushSides); }
void source_engine::bsp::File::ReadTexInfoData() { ReadData<std::remove_reference_t<decltype(m_texInfo.front())>>(LumpId::TexInfo, m_texInfo); }
void source_engine::bsp::File::ReadTexData() { ReadData<std::remove_reference_t<decltype(m_texData.front())>>(LumpId::TexData, m_texData); }
void source_engine::bsp::File::ReadModelData() { ReadData<std::remove_reference_t<decltype(m_models.front())>>(LumpId::Models, m_models); }
void source_engine::bsp::File::ReadDispLightmapSamplePositions() { ReadData<std::remove_reference_t<decltype(m_dispLightmapSamplePositions.front())>>(LumpId::DispLightmapSamplePositions, m_dispLightmapSamplePositions); }
void source_engine::bsp::File::ReadCubemapSamples() { ReadData<std::remove_reference_t<decltype(m_cubemapSamples.front())>>(LumpId::Cubemaps, m_cubemapSamples); }
void source_engine::bsp::File::ReadTexDataStringTableData() { ReadData<std::remove_reference_t<decltype(m_texDataStringTableData.front())>>(LumpId::TexDataStringTable, m_texDataStringTableData); }
void source_engine::bsp::File::ReadTexDataStringData()
{
	const auto lumpId = LumpId::TexDataStringData;
	if(HasReadLump(lumpId) == true)
		return;
	MarkLumpAsRead(lumpId);
	auto &f = m_file;
	auto &header = m_header;
	auto &lump = header.lumps.at(pragma::math::to_integral(lumpId));
	if(lump.fileLength == 0)
		return;
	f->Seek(lump.fileOffset);

	size_t offset;
	size_t size;
	auto fl = get_lump_file(lump, f, offset, size);
	if(!fl)
		return;

	while(fl->Eof() == false && fl->Tell() < offset + size) {
		m_texDataStringDataIndexMap.insert(std::make_pair(fl->Tell() - offset, m_texDataStringData.size()));
		m_texDataStringData.push_back(pragma::util::Path::CreateFile(fl->ReadString()).GetString());
	}
}
void source_engine::bsp::File::ReadNodes() { ReadData<std::remove_reference_t<decltype(m_nodes.front())>>(LumpId::Nodes, m_nodes); }
void source_engine::bsp::File::ReadLeaves()
{
	auto padding = 0u;
	if(m_header.version <= 19)
		padding = 56u - sizeof(decltype(m_leaves.front()));
	ReadData<std::remove_reference_t<decltype(m_leaves.front())>>(LumpId::Leaves, m_leaves, padding);
}
void source_engine::bsp::File::ReadLeafFaces() { ReadData<std::remove_reference_t<decltype(m_leafFaces.front())>>(LumpId::LeafFaces, m_leafFaces); }
void source_engine::bsp::File::ReadLeafBrushes() { ReadData<std::remove_reference_t<decltype(m_leafBrushes.front())>>(LumpId::LeafBrushes, m_leafBrushes); }

void source_engine::bsp::File::ReadPakfile()
{
	const auto lumpId = LumpId::PakFile;
	if(HasReadLump(lumpId) == true)
		return;
	MarkLumpAsRead(lumpId);
	auto &f = m_file;
	auto &header = m_header;
	auto &lump = header.lumps.at(pragma::math::to_integral(lumpId));
	if(lump.fileLength == 0)
		return;

	int64_t offset = lump.fileOffset + lump.fileLength - sizeof(ZipEndOfCentralDirRecord);

	auto &data = m_pakZipData;
	m_file->Seek(lump.fileOffset);
	data.resize(lump.fileLength);
	f->Read(data.data(), data.size());
	std::string err;
	m_pakZipFile = uzip::ZIPFile::Open(data.data(), data.size(), err);

	auto bFoundRecord = false;
	for(; offset >= 0; offset--) {
		f->Seek(offset);
		f->Read(&m_zipDirRecord, sizeof(m_zipDirRecord));
		constexpr uint32_t expectedSignature = (((6) << 24) | ((5) << 16) | ('K' << 8) | 'P');
		if(m_zipDirRecord.signature == expectedSignature) {
			bFoundRecord = true;
			break;
		}
	}
	if(bFoundRecord == false)
		return;
	f->Seek(lump.fileOffset + m_zipDirRecord.startOfCentralDirOffset);
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

		f->Seek(lump.fileOffset + hd.relativeOffsetOfLocalHeader);
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
		if(face.dispInfo == -1)
			continue;
		auto it = std::find_if(displacements.begin(), displacements.end(), [&face](const Displacement &disp) { return (disp.face.dispInfo == face.dispInfo) ? true : false; });
		if(it != displacements.end())
			continue;
		auto &plane = planes.at(face.planeId);
		auto &faceDispInfo = dispInfo.at(face.dispInfo);
		if(displacements.size() == displacements.capacity())
			displacements.reserve(displacements.size() + 100);
		displacements.push_back({faces.at(faceDispInfo.mapFace), plane, faceDispInfo});
	}

	auto &header = m_header;
	size_t offsetV;
	size_t sizeV;
	auto &lumpV = header.lumps.at(pragma::math::to_integral(LumpId::DispVerts));
	m_file->Seek(lumpV.fileOffset);
	auto flV = get_lump_file(lumpV, m_file, offsetV, sizeV);

	size_t offsetT;
	size_t sizeT;
	auto &lumpT = header.lumps.at(pragma::math::to_integral(LumpId::DispTris));
	m_file->Seek(lumpT.fileOffset);
	auto flT = get_lump_file(lumpT, m_file, offsetT, sizeT);
	if(!flV || !flT)
		return;
	for(auto &disp : displacements) {
		auto &dispInfo = disp.dispInfo;
		auto &dispVerts = disp.verts;
		auto &dispTris = disp.tris;
		auto p = static_cast<uint32_t>(dispInfo.power);
		auto numVerts = pragma::math::pow2(pragma::math::pow(2u, p) + 1u);
		auto numTris = 2 * pragma::math::pow2(pragma::math::pow(2u, p));

		if(sizeV != 0) {
			dispVerts.resize(numVerts);
			flV->Seek(offsetV + dispInfo.dispVertStart * sizeof(dispVerts.front()));
			flV->Read(dispVerts.data(), dispVerts.size() * sizeof(dispVerts.front()));
		}

		if(sizeT != 0) {
			dispTris.resize(numTris);
			flT->Seek(offsetT + dispInfo.dispTriStart * sizeof(dispTris.front()));
			flT->Read(dispTris.data(), dispTris.size() * sizeof(dispTris.front()));
		}
	}
}
void source_engine::bsp::File::ReadDispInfo() { ReadData<std::remove_reference_t<decltype(m_dispInfo.front())>>(LumpId::DispInfo, m_dispInfo); }

const source_engine::bsp::Lump *source_engine::bsp::File::GetLumpHeaderInfo(LumpId lumpId) const
{
	auto &header = m_header;
	return (pragma::math::to_integral(lumpId) < header.lumps.size()) ? &header.lumps.at(pragma::math::to_integral(lumpId)) : nullptr;
}
const std::vector<source_engine::bsp::GameLump> &source_engine::bsp::File::GetGameLumps()
{
	ReadGameData();
	return m_gameLumps;
}
const std::vector<source_engine::bsp::EntityBlock> &source_engine::bsp::File::GetEntities()
{
	ReadEntityData();
	return m_entities;
}
const std::vector<source_engine::bsp::Plane> &source_engine::bsp::File::GetPlanes()
{
	ReadPlaneData();
	return m_planes;
}
const std::vector<Vector3> &source_engine::bsp::File::GetVertices()
{
	ReadVertexData();
	return m_vertices;
}
const std::vector<source_engine::bsp::Edge> &source_engine::bsp::File::GetEdges()
{
	ReadEdgeData();
	return m_edges;
}
const std::vector<int32_t> &source_engine::bsp::File::GetSurfEdges()
{
	ReadSurfEdgeData();
	return m_surfEdges;
}
const std::vector<source_engine::bsp::Face> &source_engine::bsp::File::GetFaces()
{
	ReadFaceData();
	return m_faces;
}
const std::vector<source_engine::bsp::Face> &source_engine::bsp::File::GetHDRFaces()
{
	ReadHDRFaceData();
	return m_hdrFaces;
}
const std::vector<source_engine::bsp::Face> &source_engine::bsp::File::GetOriginalFaces()
{
	ReadOriginalFaceData();
	return m_origFaces;
}
const std::vector<source_engine::bsp::Brush> &source_engine::bsp::File::GetBrushes()
{
	ReadBrushData();
	return m_brushes;
}
const std::vector<source_engine::bsp::BrushSide> &source_engine::bsp::File::GetBrushSides()
{
	ReadBrushSideData();
	return m_brushSides;
}
const std::vector<source_engine::bsp::TexInfo> &source_engine::bsp::File::GetTexInfo()
{
	ReadTexInfoData();
	return m_texInfo;
}
const std::vector<source_engine::bsp::TexData> &source_engine::bsp::File::GetTexData()
{
	ReadTexData();
	return m_texData;
}
const std::vector<source_engine::bsp::Model> &source_engine::bsp::File::GetModels()
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
		if(pragma::string::compare(relFname.c_str(), "materials/", false, 10))
			relFname = pragma::string::substr(relFname, 10);

		auto it = std::find_if(strings.begin(), strings.end(), [&relFname](const std::string &fnameOther) { return pragma::string::compare(relFname, fnameOther, false); });
		if(it == strings.end())
			continue;
		std::string ext;
		if(ufile::get_extension(fname, &ext) == true && pragma::string::compare<std::string>(ext, "vmt", false) == true) {
			std::vector<uint8_t> data;
			if(ReadFile(fname, data) == true) {
				VTFLib::CVMTFile vmt {};
				if(vmt.Load(data.data(), static_cast<vlUInt>(data.size())) == vlTrue) {
					auto *vmtRoot = vmt.GetRoot();
					std::string shader = vmtRoot->GetName();
					if(pragma::string::compare<std::string>(shader, "patch", false)) {
						VTFLib::Nodes::CVMTNode *node = nullptr;
						if((node = vmtRoot->GetNode("include")) != nullptr) {
							if(node->GetType() == VMTNodeType::NODE_TYPE_STRING) {
								auto *includeNode = static_cast<VTFLib::Nodes::CVMTStringNode *>(node);
								std::string include = includeNode->GetValue();
								ufile::remove_extension_from_filename(include);
								if(pragma::string::compare(include.c_str(), "materials/", false, 10))
									*it = pragma::string::substr(include, 10);
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
const std::vector<source_engine::bsp::Node> &source_engine::bsp::File::GetNodes()
{
	ReadNodes();
	return m_nodes;
}
const std::vector<source_engine::bsp::Leaf> &source_engine::bsp::File::GetLeaves()
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
const std::vector<source_engine::bsp::DisplacementInfo> &source_engine::bsp::File::GetDispInfo()
{
	ReadDispInfo();
	return m_dispInfo;
}
const std::vector<source_engine::bsp::Displacement> &source_engine::bsp::File::GetDisplacements()
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
const std::vector<source_engine::bsp::CubemapSample> &source_engine::bsp::File::GetCubemapSamples()
{
	ReadCubemapSamples();
	return m_cubemapSamples;
}
const source_engine::bsp::StaticPropData &source_engine::bsp::File::GetStaticPropData()
{
	ReadStaticPropsData();
	return m_staticPropData;
}
const source_engine::bsp::Header &source_engine::bsp::File::GetHeaderData() const { return m_header; }
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
		return pragma::string::compare(fname,fnameOther,false);
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
