// SPDX-FileCopyrightText: (c) 2024 Silverlan <opensource@pragma-engine.com>
// SPDX-License-Identifier: MIT

export module source_engine.bsp;

import util_zip;
export import source_engine.vmf;

export namespace source_engine::bsp {
	constexpr uint8_t HEADER_LUMPS = 64u;
	constexpr uint8_t MAX_DISP_CORNER_NEIGHBORS = 4;
	enum class LumpId : uint8_t {
		Entities = 0u,
		Planes = 1u,
		TexData = 2u,
		Vertices = 3u,
		Visibility = 4u,
		Nodes = 5u,
		TexInfo = 6u,
		Faces = 7u,
		Lighting = 8u,
		Leaves = 10u,
		Edges = 12u,
		SurfEdges = 13u,
		Models = 14u,
		LeafFaces = 16u,
		LeafBrushes = 17u,
		Brushes = 18u,
		BrushSides = 19u,
		DispInfo = 26u,
		OriginalFaces = 27u,
		DispVerts = 33u,
		DispLightmapSamplePositions = 34u,
		Game = 35u,
		PakFile = 40u,
		Cubemaps = 42u,
		TexDataStringData = 43u,
		TexDataStringTable = 44u,
		DispTris = 48u,
		LightingHdr = 53u,
		FacesHdr = 58u,
	};

	struct Leaf {
		int32_t contents;
		int16_t cluster;
		int16_t area : 9;
		int16_t flags : 7;
		std::array<int16_t, 3> mins;
		std::array<int16_t, 3> maxs;
		uint16_t firstLeafFace;
		uint16_t numLeafFaces;
		uint16_t firstLeafBrush;
		uint16_t numLeafBrushes;
		int16_t leafWaterDataId;
	};
#pragma pack(push, 1)
	struct ColorRgbExp32 {
		uint8_t r, g, b;
		int8_t exponent;
	};
	struct Lump {
		int32_t fileOffset;
		int32_t fileLength;
		int32_t version;
		std::array<int8_t, 4> fourCc;
	};
	struct Header {
		int32_t identifier;
		int32_t version;
		std::array<Lump, HEADER_LUMPS> lumps;
		int32_t mapRevision;
	};
	struct Plane {
		Vector3 normal;
		float dist;
		int32_t type;
	};
	struct Edge {
		std::array<uint16_t, 2> vertexIndices;
	};
	struct Face {
		uint16_t planeId;
		uint8_t side;
		uint8_t onNode;
		int32_t firstEdge;
		int16_t numEdges;
		int16_t texInfo;
		int16_t dispInfo;
		int16_t surfaceFogVolumeId;
		std::array<uint8_t, 4> styles;
		int32_t lightmapLumpOffset;
		float area;
		std::array<int32_t, 2> lightmapTextureMinsInLuxels;
		std::array<int32_t, 2> lightmapTextureSizeInLuxels;
		int32_t origFace;
		uint16_t numPrimitives;
		uint16_t firstPrimId;
		uint32_t smoothingGroups;
	};
	struct Brush {
		int32_t firstSide;
		int32_t numSides;
		int32_t contentFlags;
	};
	struct BrushSide {
		uint16_t planeId;
		int16_t texInfo;
		int16_t dispInfo;
		int16_t bevel;
	};
	struct TexInfo {
		std::array<std::array<float, 4>, 2> textureVecs;
		std::array<std::array<float, 4>, 2> lightmapVecs;
		int32_t flags;
		int32_t texData;
	};
	struct TexData {
		Vector3 reflectivity;
		int32_t nameStringTableId;
		int32_t width, height;
		int32_t viewWidth, viewHeight;
	};
	struct Model {
		Vector3 mins, maxs;
		Vector3 origin;
		int32_t headNode;
		int32_t firstFace, numFaces;
	};
	struct Node {
		int32_t planeIndex;
		std::array<int32_t, 2> children;
		std::array<int16_t, 3> mins;
		std::array<int16_t, 3> maxs;
		uint16_t firstFace;
		uint16_t numFaces;
		int16_t area;
		int16_t paddding;
	};
	struct CubemapSample {
		std::array<int32_t, 3> origin;
		int32_t resolution;
	};
#pragma pack(pop)
	enum class NeighborSpan : uint8_t { CornerToCorner = 0, CornerToMidpoint, MidpointToCorner };
	enum class NeighborOrientation : uint8_t {
		OrientationCcw0 = 0,
		OrientationCcw90,
		OrientationCcw180,
		OrientationCcw270,
	};
	struct DisplacementSubNeighbor {
	  public:
		uint16_t GetNeighborIndex() const { return m_neighborDispInfoIndex; }
		NeighborSpan GetSpan() const { return m_span; }
		NeighborSpan GetNeighborSpan() const { return m_neighborSpan; }
		NeighborOrientation GetNeighborOrientation() const { return m_neighborOrientation; }

		bool IsValid() const { return m_neighborDispInfoIndex != 0xFFFF; }
		void SetInvalid() { m_neighborDispInfoIndex = 0xFFFF; }
	  public:
		uint16_t m_neighborDispInfoIndex;

		NeighborOrientation m_neighborOrientation;

		NeighborSpan m_span;
		NeighborSpan m_neighborSpan;
	};
#pragma pack(push, 1)
	class DisplacementNeighbor {
	  public:
		void SetInvalid()
		{
			m_subNeighbors[0].SetInvalid();
			m_subNeighbors[1].SetInvalid();
		}

		bool IsValid() { return m_subNeighbors[0].IsValid() || m_subNeighbors[1].IsValid(); }
	  public:
		std::array<DisplacementSubNeighbor, 2> m_subNeighbors;
	};
#pragma pack(pop)
	class DisplacementCornerNeighbors {
	  public:
		void SetInvalid() { m_nNeighbors = 0; }
	  public:
		std::array<uint16_t, MAX_DISP_CORNER_NEIGHBORS> m_Neighbors;
		uint8_t m_nNeighbors;
	};
	struct DisplacementInfo {
		Vector3 startPosition;
		int32_t dispVertStart;
		int32_t dispTriStart;
		int32_t power;
		int32_t minTess;
		float smoothingAngle;
		int32_t contents;
		uint16_t mapFace;
		int32_t lightmapAlphaStart;
		int32_t lightmapSamplePositionStart;
		std::array<DisplacementNeighbor, 4> edgeNeighbors;
		std::array<DisplacementCornerNeighbors, 4> cornerNeighbors;
		std::array<uint32_t, 10> allowedVerts;
	};
#pragma pack(push, 1)
	struct DisplacementVertex {
		Vector3 vec;
		float dist;
		float alpha;
	};
	struct DisplacementTriangle {
		uint16_t Tags;
	};
	struct ZipFileHeader {
		uint32_t signature;
		uint16_t versionMadeBy;
		uint16_t versionNeededToExtract;
		uint16_t flags;
		uint16_t compressionMethod;
		uint16_t lastModifiedTime;
		uint16_t lastModifiedDate;
		uint32_t crc32;
		uint32_t compressedSize;
		uint32_t uncompressedSize;
		uint16_t fileNameLength;
		uint16_t extraFieldLength;
		uint16_t fileCommentLength;
		uint16_t diskNumberStart;
		uint16_t internalFileAttribs;
		uint32_t externalFileAttribs;
		uint32_t relativeOffsetOfLocalHeader;
	};
	struct ZipLocalFileHeader {
		uint32_t signature;
		uint16_t versionNeededToExtract;
		uint16_t flags;
		uint16_t compressionMethod;
		uint16_t lastModifiedTime;
		uint16_t lastModifiedDate;
		uint32_t crc32;
		uint32_t compressedSize;
		uint32_t uncompressedSize;
		uint16_t fileNameLength;
		uint16_t extraFieldLength;
	};
	struct ZipEndOfCentralDirRecord {
		uint32_t signature;
		uint16_t numberOfThisDisk;
		uint16_t numberOfTheDiskWithStartOfCentralDirectory;
		uint16_t nCentralDirectoryEntries_ThisDisk;
		uint16_t nCentralDirectoryEntries_Total;
		uint32_t centralDirectorySize;
		uint32_t startOfCentralDirOffset;
		uint16_t commentLength;
	};
	struct GameLump {
		int32_t id;
		uint16_t flags;
		uint16_t version;
		int32_t fileOffset;
		int32_t fileLength;
	};
#pragma pack(pop)
	struct StaticPropLump {
		Vector3 origin;
		EulerAngles angles;
		uint16_t propType;
		uint16_t firstLeaf;
		uint16_t leafCount;
		uint8_t solid;
		uint8_t flags;
		int32_t skin;
		float fadeMinDist;
		float fadeMaxDist;
		Vector3 lightingOrigin;
	};
	struct StaticPropData {
		std::vector<std::string> dictionaryModelNames {};
		std::vector<uint16_t> leaves {};
		std::vector<StaticPropLump> staticPropLumps {};
	};
	struct Displacement {
		Displacement(const Face &_face, const Plane &_plane, const DisplacementInfo &_dispInfo) : face(_face), plane(_plane), dispInfo(_dispInfo) {}
		std::vector<DisplacementVertex> verts;
		std::vector<DisplacementTriangle> tris;
		const Face &face;
		const Plane &plane;
		const DisplacementInfo &dispInfo;
	};

	enum class ResultCode : uint32_t { Success = 0, FileNotFound, InvalidHeaderIdent };

	using EntityBlock = std::shared_ptr<vmf::DataFileBlock>;
	class File {
	  public:
		static std::unique_ptr<File> Open(pragma::fs::VFilePtr &f, ResultCode &code);

		const Lump *GetLumpHeaderInfo(LumpId lumpId) const;
		const std::vector<GameLump> &GetGameLumps();
		const std::vector<EntityBlock> &GetEntities();
		const std::vector<Plane> &GetPlanes();
		const std::vector<Vector3> &GetVertices();
		const std::vector<Edge> &GetEdges();
		const std::vector<int32_t> &GetSurfEdges();
		const std::vector<Face> &GetFaces();
		const std::vector<Face> &GetHDRFaces();
		const std::vector<Face> &GetOriginalFaces();
		const std::vector<Brush> &GetBrushes();
		const std::vector<BrushSide> &GetBrushSides();
		const std::vector<TexInfo> &GetTexInfo();
		const std::vector<TexData> &GetTexData();
		const std::vector<Model> &GetModels();
		const std::vector<uint32_t> &GetTexDataStringIndices();
		const std::vector<std::string> &GetTexDataStrings();
		const std::vector<std::string> &GetTranslatedTexDataStrings();
		const std::vector<Node> &GetNodes();
		const std::vector<Leaf> &GetLeaves();
		const std::vector<uint16_t> &GetLeafFaces();
		const std::vector<uint16_t> &GetLeafBrushes();
		const std::vector<std::string> &GetFilenames();
		const std::vector<DisplacementInfo> &GetDispInfo();
		const std::vector<Displacement> &GetDisplacements();
		const std::vector<uint8_t> &GetLightMapData();
		const std::vector<uint8_t> &GetHDRLightMapData();
		const std::vector<std::vector<uint8_t>> &GetVisibilityData();
		const std::vector<uint8_t> &GetDispLightmapSamplePositions();
		const std::vector<CubemapSample> &GetCubemapSamples();
		const StaticPropData &GetStaticPropData();
		const Header &GetHeaderData() const;
		const bool ReadFile(const std::string &fname, std::vector<uint8_t> &data);
	  private:
		pragma::fs::VFilePtr m_file;
		Header m_header;
		uint64_t m_readLumps = 0;
		bool m_bStringDataIndicesRead = false;
		bool m_bStringDataTranslated = false;
		bool m_bDisplacementsRead = false;
		File(pragma::fs::VFilePtr &f, const Header &header);

		std::vector<GameLump> m_gameLumps;
		std::vector<EntityBlock> m_entities;
		std::vector<Plane> m_planes;
		std::vector<Vector3> m_vertices;
		std::vector<Edge> m_edges;
		std::vector<int32_t> m_surfEdges;
		std::vector<Face> m_faces;
		std::vector<Face> m_hdrFaces;
		std::vector<Face> m_origFaces;
		std::vector<uint8_t> m_lightMapData;
		std::vector<uint8_t> m_lightMapDataHDR;
		std::vector<Brush> m_brushes;
		std::vector<BrushSide> m_brushSides;
		std::vector<TexInfo> m_texInfo;
		std::vector<TexData> m_texData;
		std::vector<Model> m_models;
		std::vector<int32_t> m_texDataStringTableData;
		std::vector<std::string> m_texDataStringData;
		std::vector<std::string> m_texDataStringDataTranslated;
		std::unordered_map<uint64_t, uint32_t> m_texDataStringDataIndexMap;
		std::vector<uint32_t> m_texDataStringDataIndices;
		std::vector<Node> m_nodes;
		std::vector<Leaf> m_leaves;
		std::vector<uint16_t> m_leafFaces;
		std::vector<uint16_t> m_leafBrushes;
		std::vector<DisplacementInfo> m_dispInfo;
		std::vector<Displacement> m_displacements;
		std::vector<std::vector<uint8_t>> m_visibilityData;
		std::vector<uint8_t> m_dispLightmapSamplePositions;
		std::vector<CubemapSample> m_cubemapSamples;
		StaticPropData m_staticPropData = {};

		std::vector<uint8_t> m_pakZipData;
		std::shared_ptr<uzip::ZIPFile> m_pakZipFile = nullptr;
		ZipEndOfCentralDirRecord m_zipDirRecord;
		std::vector<ZipFileHeader> m_fileHeaders;
		std::vector<ZipLocalFileHeader> m_localFileHeaders;
		std::vector<uint64_t> m_fileDataOffsets;
		std::vector<std::string> m_fileNames;

		bool HasReadLump(LumpId lumpId) const;
		void MarkLumpAsRead(LumpId lumpId);

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
		void ReadData(LumpId lumpId, TContainer &data);
		template<class T, class TContainer>
		void ReadData(LumpId lumpId, TContainer &data, uint32_t padding);
	};
};

namespace source_engine::bsp {
	bool lzma_uncompress(unsigned char *dest, size_t *destLen, const unsigned char *src, size_t *srcLen, const unsigned char *props, size_t propsSize);
};
