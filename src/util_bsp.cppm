// SPDX-FileCopyrightText: (c) 2024 Silverlan <opensource@pragma-engine.com>
// SPDX-License-Identifier: MIT

module;

#include <cinttypes>
#include <array>
#include <unordered_map>
#include <mathutil/uvec.h>
#include <fsys/vfileptr.h>

#define MAX_DISP_CORNER_NEIGHBORS 4

export module source_engine.bsp;

import util_zip;
export import source_engine.vmf;

export namespace source_engine::bsp {
	const auto HEADER_LUMPS = 64u;

	const auto LUMP_ID_ENTITIES = 0u;
	const auto LUMP_ID_PLANES = 1u;
	const auto LUMP_ID_TEXDATA = 2u;
	const auto LUMP_ID_VERTICES = 3u;
	const auto LUMP_ID_VISIBILITY = 4u;
	const auto LUMP_ID_TEXINFO = 6u;
	const auto LUMP_ID_FACES = 7u;
	const auto LUMP_ID_LIGHTING = 8u;
	const auto LUMP_ID_EDGES = 12u;
	const auto LUMP_ID_SURF_EDGES = 13u;
	const auto LUMP_ID_MODELS = 14u;
	const auto LUMP_ID_BRUSHES = 18u;
	const auto LUMP_ID_BRUSH_SIDES = 19u;
	const auto LUMP_ID_ORIGINAL_FACES = 27u;
	const auto LUMP_ID_DISP_VERTS = 33u;
	const auto LUMP_ID_DISP_LIGHTMAP_SAMPLE_POSITIONS = 34u;
	const auto LUMP_ID_GAME = 35u;
	const auto LUMP_ID_PAKFILE = 40u;
	const auto LUMP_ID_CUBEMAPS = 42u;
	const auto LUMP_ID_TEXDATA_STRING_DATA = 43u;
	const auto LUMP_ID_TEXDATA_STRING_TABLE = 44u;
	const auto LUMP_ID_DISP_TRIS = 48u;
	const auto LUMP_ID_LIGHTING_HDR = 53u;
	const auto LUMP_ID_FACES_HDR = 58u;
	struct dleaf_t {
		int32_t contents;            // OR of all brushes (not needed?)
		int16_t cluster;             // cluster this leaf is in
		int16_t area : 9;            // area this leaf is in
		int16_t flags : 7;           // flags
		std::array<int16_t, 3> mins; // for frustum culling
		std::array<int16_t, 3> maxs;
		uint16_t firstleafface; // index into leaffaces
		uint16_t numleaffaces;
		uint16_t firstleafbrush; // index into leafbrushes
		uint16_t numleafbrushes;
		int16_t leafWaterDataID; // -1 for not in water

		//!!! NOTE: for maps of version 19 or lower uncomment this block
		/*
		CompressedLightCube	ambientLighting;	// Precaculated light info for entities.
		short			padding;		// padding to 4-byte boundary
		*/
	};
#pragma pack(push, 1)
	struct ColorRGBExp32 {
		uint8_t r, g, b;
		int8_t exponent;
	};
	struct lump_t {
		int32_t fileofs;              // offset into file (bytes)
		int32_t filelen;              // length of lump (bytes)
		int32_t version;              // lump format version
		std::array<int8_t, 4> fourCC; // lump ident code
	};
	struct dheader_t {
		int32_t ident;                          // BSP file identifier
		int32_t version;                        // BSP file version
		std::array<lump_t, HEADER_LUMPS> lumps; // lump directory array
		int32_t mapRevision;                    // the map's revision (iteration, version) number
	};
	struct dplane_t {
		Vector3 normal; // normal vector
		float dist;     // distance from origin
		int32_t type;   // plane axis identifier
	};
	struct dedge_t {
		std::array<uint16_t, 2> v; // vertex indices
	};
	struct dface_t {
		uint16_t planenum;                                  // the plane number
		uint8_t side;                                       // faces opposite to the node's plane direction
		uint8_t onNode;                                     // 1 of on node, 0 if in leaf
		int32_t firstedge;                                  // index into surfedges
		int16_t numedges;                                   // number of surfedges
		int16_t texinfo;                                    // texture info
		int16_t dispinfo;                                   // displacement info
		int16_t surfaceFogVolumeID;                         // ?
		std::array<uint8_t, 4> styles;                      // switchable lighting info
		int32_t lightofs;                                   // offset into lightmap lump
		float area;                                         // face area in units^2
		std::array<int32_t, 2> LightmapTextureMinsInLuxels; // texture lighting info
		std::array<int32_t, 2> LightmapTextureSizeInLuxels; // texture lighting info
		int32_t origFace;                                   // original face this was split from
		uint16_t numPrims;                                  // primitives
		uint16_t firstPrimID;
		uint32_t smoothingGroups; // lightmap smoothing group
	};
	struct dbrush_t {
		int32_t firstside; // first brushside
		int32_t numsides;  // number of brushsides
		int32_t contents;  // contents flags
	};
	struct dbrushside_t {
		uint16_t planenum; // facing out of the leaf
		int16_t texinfo;   // texture info
		int16_t dispinfo;  // displacement info
		int16_t bevel;     // is the side a bevel plane?
	};
	struct texinfo_t {
		std::array<std::array<float, 4>, 2> textureVecs;  // [s/t][xyz offset]
		std::array<std::array<float, 4>, 2> lightmapVecs; // [s/t][xyz offset] - length is in units of texels/area
		int32_t flags;                                    // miptex flags	overrides
		int32_t texdata;                                  // Pointer to texture name, size, etc.
	};
	struct dtexdata_t {
		Vector3 reflectivity;      // RGB reflectivity
		int32_t nameStringTableID; // index into TexdataStringTable
		int32_t width, height;     // source image
		int32_t view_width, view_height;
	};
	struct dmodel_t {
		Vector3 mins, maxs;          // bounding box
		Vector3 origin;              // for sounds or lights
		int32_t headnode;            // index into node array
		int32_t firstface, numfaces; // index into face array
	};
	struct dnode_t {
		int32_t planenum;                // index into plane array
		std::array<int32_t, 2> children; // negative numbers are -(leafs + 1), not nodes
		std::array<int16_t, 3> mins;     // for frustum culling
		std::array<int16_t, 3> maxs;
		uint16_t firstface; // index into face array
		uint16_t numfaces;  // counting both sides
		int16_t area;       // If all leaves below this node are in the same area, then
		// this is the area index. If not, this is -1.
		int16_t paddding; // pad to 32 bytes length
	};
	struct dcubemapsample_t {
		std::array<int32_t, 3> origin; // position of light snapped to the nearest integer
		int32_t size;                  // resolution of cubemap, 0 - default
	};
#pragma pack(pop)
	typedef enum { CORNER_TO_CORNER = 0, CORNER_TO_MIDPOINT = 1, MIDPOINT_TO_CORNER = 2 } NeighborSpan;
	typedef enum { ORIENTATION_CCW_0 = 0, ORIENTATION_CCW_90 = 1, ORIENTATION_CCW_180 = 2, ORIENTATION_CCW_270 = 3 } NeighborOrientation;
	struct CDispSubNeighbor {
	  public:
		unsigned short GetNeighborIndex() const { return m_iNeighbor; }
		NeighborSpan GetSpan() const { return (NeighborSpan)m_Span; }
		NeighborSpan GetNeighborSpan() const { return (NeighborSpan)m_NeighborSpan; }
		NeighborOrientation GetNeighborOrientation() const { return (NeighborOrientation)m_NeighborOrientation; }

		bool IsValid() const { return m_iNeighbor != 0xFFFF; }
		void SetInvalid() { m_iNeighbor = 0xFFFF; }
	  public:
		uint16_t m_iNeighbor; // This indexes into ddispinfos.
		// 0xFFFF if there is no neighbor here.

		uint8_t m_NeighborOrientation; // (CCW) rotation of the neighbor wrt this displacement.

		// These use the NeighborSpan type.
		uint8_t m_Span;         // Where the neighbor fits onto this side of our displacement.
		uint8_t m_NeighborSpan; // Where we fit onto our neighbor.
	};
#pragma pack(push, 1)
	class CDispNeighbor {
	  public:
		void SetInvalid()
		{
			m_SubNeighbors[0].SetInvalid();
			m_SubNeighbors[1].SetInvalid();
		}

		// Returns false if there isn't anything touching this edge.
		bool IsValid() { return m_SubNeighbors[0].IsValid() || m_SubNeighbors[1].IsValid(); }
	  public:
		// Note: if there is a neighbor that fills the whole side (CORNER_TO_CORNER),
		//       then it will always be in CDispNeighbor::m_Neighbors[0]
		std::array<CDispSubNeighbor, 2> m_SubNeighbors;
	};
#pragma pack(pop)
	class CDispCornerNeighbors {
	  public:
		void SetInvalid() { m_nNeighbors = 0; }
	  public:
		std::array<uint16_t, MAX_DISP_CORNER_NEIGHBORS> m_Neighbors; // indices of neighbors.
		uint8_t m_nNeighbors;
	};
	struct ddispinfo_t {
		Vector3 startPosition;                               // start position used for orientation
		int32_t DispVertStart;                               // Index into LUMP_DISP_VERTS.
		int32_t DispTriStart;                                // Index into LUMP_DISP_TRIS.
		int32_t power;                                       // power - indicates size of surface (2^power	1)
		int32_t minTess;                                     // minimum tesselation allowed
		float smoothingAngle;                                // lighting smoothing angle
		int32_t contents;                                    // surface contents
		uint16_t MapFace;                                    // Which map face this displacement comes from.
		int32_t LightmapAlphaStart;                          // Index into ddisplightmapalpha.
		int32_t LightmapSamplePositionStart;                 // Index into LUMP_DISP_LIGHTMAP_SAMPLE_POSITIONS.
		std::array<CDispNeighbor, 4> EdgeNeighbors;          // Indexed by NEIGHBOREDGE_ defines.
		std::array<CDispCornerNeighbors, 4> CornerNeighbors; // Indexed by CORNER_ defines.
		std::array<uint32_t, 10> AllowedVerts;               // active verticies
	};
#pragma pack(push, 1)
	struct dDispVert {
		Vector3 vec; // Vector field defining displacement volume.
		float dist;  // Displacement distances.
		float alpha; // "per vertex" alpha values.
	};
	struct dDispTri {
		uint16_t Tags; // Displacement triangle tags.
	};
	struct ZIP_FileHeader {
		uint32_t signature;                   //  4 bytes PK12
		uint16_t versionMadeBy;               // version made by 2 bytes
		uint16_t versionNeededToExtract;      // version needed to extract 2 bytes
		uint16_t flags;                       // general purpose bit flag 2 bytes
		uint16_t compressionMethod;           // compression method 2 bytes
		uint16_t lastModifiedTime;            // last mod file time 2 bytes
		uint16_t lastModifiedDate;            // last mod file date 2 bytes
		uint32_t crc32;                       // crc-32 4 bytes
		uint32_t compressedSize;              // compressed size 4 bytes
		uint32_t uncompressedSize;            // uncompressed size 4 bytes
		uint16_t fileNameLength;              // file name length 2 bytes
		uint16_t extraFieldLength;            // extra field length 2 bytes
		uint16_t fileCommentLength;           // file comment length 2 bytes
		uint16_t diskNumberStart;             // disk number start 2 bytes
		uint16_t internalFileAttribs;         // internal file attributes 2 bytes
		uint32_t externalFileAttribs;         // external file attributes 4 bytes
		uint32_t relativeOffsetOfLocalHeader; // relative offset of local header 4 bytes
		                                      // file name (variable size)
		                                      // extra field (variable size)
		                                      // file comment (variable size)
	};
	struct ZIP_LocalFileHeader {
		uint32_t signature;              //local file header signature 4 bytes PK34
		uint16_t versionNeededToExtract; // version needed to extract 2 bytes
		uint16_t flags;                  // general purpose bit flag 2 bytes
		uint16_t compressionMethod;      // compression method 2 bytes
		uint16_t lastModifiedTime;       // last mod file time 2 bytes
		uint16_t lastModifiedDate;       // last mod file date 2 bytes
		uint32_t crc32;                  // crc-32 4 bytes
		uint32_t compressedSize;         // compressed size 4 bytes
		uint32_t uncompressedSize;       // uncompressed size 4 bytes
		uint16_t fileNameLength;         // file name length 2 bytes
		uint16_t extraFieldLength;       // extra field length 2 bytes
		                                 // file name (variable size)
		                                 // extra field (variable size)
	};
	struct ZIP_EndOfCentralDirRecord {
		uint32_t signature;                                  // 4 bytes PK56
		uint16_t numberOfThisDisk;                           // 2 bytes
		uint16_t numberOfTheDiskWithStartOfCentralDirectory; // 2 bytes
		uint16_t nCentralDirectoryEntries_ThisDisk;          // 2 bytes
		uint16_t nCentralDirectoryEntries_Total;             // 2 bytes
		uint32_t centralDirectorySize;                       // 4 bytes
		uint32_t startOfCentralDirOffset;                    // 4 bytes
		uint16_t commentLength;                              // 2 bytes
		                                                     // zip file comment follows
	};
	struct dgamelump_t {
		int32_t id;       // gamelump ID
		uint16_t flags;   // flags
		uint16_t version; // gamelump version
		int32_t fileofs;  // offset to this gamelump
		int32_t filelen;  // length
	};
#pragma pack(pop)
	struct StaticPropLump_t {
		Vector3 Origin;     // origin
		EulerAngles Angles; // orientation (pitch roll yaw)
		uint16_t PropType;  // index into model name dictionary
		uint16_t FirstLeaf; // index into leaf array
		uint16_t LeafCount;
		uint8_t Solid; // solidity type
		uint8_t Flags;
		int32_t Skin; // model skin numbers
		float FadeMinDist;
		float FadeMaxDist;
		Vector3 LightingOrigin; // for lighting
	};
	struct StaticPropData {
		std::vector<std::string> dictionaryModelNames {};
		std::vector<uint16_t> leaves {};
		std::vector<StaticPropLump_t> staticPropLumps {};
	};
	struct dDisp {
		dDisp(const dface_t &_face, const dplane_t &_plane, const ddispinfo_t &_dispInfo) : face(_face), plane(_plane), dispInfo(_dispInfo) {}
		std::vector<dDispVert> verts;
		std::vector<dDispTri> tris;
		const dface_t &face;
		const dplane_t &plane;
		const ddispinfo_t &dispInfo;
	};

	enum class ResultCode : uint32_t { Success = 0, FileNotFound, InvalidHeaderIdent };

	using EntityBlock = std::shared_ptr<vmf::DataFileBlock>;
	class File {
	  public:
		static std::unique_ptr<File> Open(VFilePtr &f, ResultCode &code);

		const lump_t *GetLumpHeaderInfo(uint32_t lumpId) const;
		const std::vector<dgamelump_t> &GetGameLumps();
		const std::vector<EntityBlock> &GetEntities();
		const std::vector<dplane_t> &GetPlanes();
		const std::vector<Vector3> &GetVertices();
		const std::vector<dedge_t> &GetEdges();
		const std::vector<int32_t> &GetSurfEdges();
		const std::vector<dface_t> &GetFaces();
		const std::vector<dface_t> &GetHDRFaces();
		const std::vector<dface_t> &GetOriginalFaces();
		const std::vector<dbrush_t> &GetBrushes();
		const std::vector<dbrushside_t> &GetBrushSides();
		const std::vector<texinfo_t> &GetTexInfo();
		const std::vector<dtexdata_t> &GetTexData();
		const std::vector<dmodel_t> &GetModels();
		const std::vector<uint32_t> &GetTexDataStringIndices();
		const std::vector<std::string> &GetTexDataStrings();
		const std::vector<std::string> &GetTranslatedTexDataStrings();
		const std::vector<dnode_t> &GetNodes();
		const std::vector<dleaf_t> &GetLeaves();
		const std::vector<uint16_t> &GetLeafFaces();
		const std::vector<uint16_t> &GetLeafBrushes();
		const std::vector<std::string> &GetFilenames();
		const std::vector<ddispinfo_t> &GetDispInfo();
		const std::vector<dDisp> &GetDisplacements();
		const std::vector<uint8_t> &GetLightMapData();
		const std::vector<uint8_t> &GetHDRLightMapData();
		const std::vector<std::vector<uint8_t>> &GetVisibilityData();
		const std::vector<uint8_t> &GetDispLightmapSamplePositions();
		const std::vector<dcubemapsample_t> &GetCubemapSamples();
		const StaticPropData &GetStaticPropData();
		const dheader_t &GetHeaderData() const;
		const bool ReadFile(const std::string &fname, std::vector<uint8_t> &data);
	  private:
		VFilePtr m_file;
		dheader_t m_header;
		uint64_t m_readLumps = 0;
		bool m_bStringDataIndicesRead = false;
		bool m_bStringDataTranslated = false;
		bool m_bDisplacementsRead = false;
		File(VFilePtr &f, const dheader_t &header);

		std::vector<dgamelump_t> m_gameLumps;
		std::vector<EntityBlock> m_entities;
		std::vector<dplane_t> m_planes;
		std::vector<Vector3> m_vertices;
		std::vector<dedge_t> m_edges;
		std::vector<int32_t> m_surfEdges;
		std::vector<dface_t> m_faces;
		std::vector<dface_t> m_hdrFaces;
		std::vector<dface_t> m_origFaces;
		std::vector<uint8_t> m_lightMapData;
		std::vector<uint8_t> m_lightMapDataHDR;
		std::vector<dbrush_t> m_brushes;
		std::vector<dbrushside_t> m_brushSides;
		std::vector<texinfo_t> m_texInfo;
		std::vector<dtexdata_t> m_texData;
		std::vector<dmodel_t> m_models;
		std::vector<int32_t> m_texDataStringTableData;
		std::vector<std::string> m_texDataStringData;
		std::vector<std::string> m_texDataStringDataTranslated;
		std::unordered_map<uint64_t, uint32_t> m_texDataStringDataIndexMap;
		std::vector<uint32_t> m_texDataStringDataIndices;
		std::vector<dnode_t> m_nodes;
		std::vector<dleaf_t> m_leaves;
		std::vector<uint16_t> m_leafFaces;
		std::vector<uint16_t> m_leafBrushes;
		std::vector<ddispinfo_t> m_dispInfo;
		std::vector<dDisp> m_displacements;
		std::vector<std::vector<uint8_t>> m_visibilityData;
		std::vector<uint8_t> m_dispLightmapSamplePositions;
		std::vector<dcubemapsample_t> m_cubemapSamples;
		StaticPropData m_staticPropData = {};

		std::vector<uint8_t> m_pakZipData;
		std::shared_ptr<uzip::ZIPFile> m_pakZipFile = nullptr;
		ZIP_EndOfCentralDirRecord m_zipDirRecord;
		std::vector<ZIP_FileHeader> m_fileHeaders;
		std::vector<ZIP_LocalFileHeader> m_localFileHeaders;
		std::vector<uint64_t> m_fileDataOffsets;
		std::vector<std::string> m_fileNames;

		bool HasReadLump(uint32_t lumpId) const;
		void MarkLumpAsRead(uint32_t lumpId);

		void ReadGameData();
		void ReadStaticPropsData();
		void ReadEntityData();
		void ReadPlaneData();
		void ReadVertexData();
		void ReadEdgeData();
		void ReadSurfEdgeData();
		void ReadFaceData();
		void ReadHDRFaceData();
		void ReadOriginalFaceData();
		void ReadBrushData();
		void ReadBrushSideData();
		void ReadTexInfoData();
		void ReadTexData();
		void ReadModelData();
		void ReadTexDataStringTableData();
		void ReadTexDataStringData();
		void ReadNodes();
		void ReadLeaves();
		void ReadLeafFaces();
		void ReadLeafBrushes();
		void ReadPakfile();
		void ReadDisplacementData();
		void ReadDispInfo();
		void ReadLightMapData();
		void ReadHDRLightMapData();
		void ReadVisibilityData();
		void ReadDispLightmapSamplePositions();
		void ReadCubemapSamples();
		template<class T, class TContainer>
		void ReadData(uint32_t lumpId, TContainer &data);
		template<class T, class TContainer>
		void ReadData(uint32_t lumpId, TContainer &data, uint32_t padding);
	};
};

namespace source_engine::bsp {
	bool lzma_uncompress(unsigned char *dest, size_t *destLen, const unsigned char *src, size_t *srcLen, const unsigned char *props, size_t propsSize);
};
