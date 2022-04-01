#include <iomanip>
#include <pre/geometry/Bound>
#include <pre/geometry/Triangle>
#include <pre/random>
#include <pre/terminal>
#include <regex>

#include <H5Cpp.h>

using namespace std::string_literals;
#define THROW_IF(what, message)                                              \
  do {                                                                       \
    if ((what))                                                              \
      throw std::runtime_error(std::string(message) + " ["s + #what + "]"s); \
  } while (false)

const char* programName = nullptr;

bool alwaysYes = false;

void QuitWithHint() {
  std::cerr << "Type '" << programName << " -h' for program help.\n";
  std::exit(EXIT_FAILURE);
}

struct VoxelRegion {
  enum class Mode { Area, CosineArea, Normal };
  void add(pre::Vec3f point, float area) {
    pre::Vec3<uint64_t> index = (point - bound.lower()) * invVoxelSize;
    if (!(index < count).all()) return;
    voxelData[(index[0] * count[1] + index[1]) * count[2] + index[2]] += area;
  }
  void add(pre::Vec3f point, pre::Vec3f normal) {
    pre::Vec3<uint64_t> index = (point - bound.lower()) * invVoxelSize;
    if (!(index < count).all()) return;
    uint64_t i = (index[0] * count[1] + index[1]) * count[2] + index[2];
    voxelData[i * 3 + 0] += normal[0];
    voxelData[i * 3 + 1] += normal[1];
    voxelData[i * 3 + 2] += normal[2];
  }
  void add(const pre::Triangle3f& tri) {
    pre::Bound3f triBound = tri;
    if (!bound.overlaps(triBound)) return;
    pre::Vec3f normal = tri.normal();
    if (normal[2] < 0) normal *= -1;
    float area = pre::sqrt(pre::length2(normal)) * 0.5f;
    if (mode == Mode::CosineArea) area = normal[2];
    if (!std::isfinite(area)) return;
    pre::Vec3f triExtent = triBound.extent();
    pre::Vec3f triCenter = tri.center();
    if ((triExtent < 0.25f * voxelSize).all()) {
      if (mode == Mode::Normal)
        add(triCenter, normal);
      else
        add(triCenter, area);
    } else {
      pre::Vec3f ratio = triExtent * invVoxelSize;
      int oversample = 32 * pre::max(1.0f, ratio[ratio.argmax()]);
      normal /= oversample;
      area /= oversample;
      while (oversample-- > 0) {
        auto [u0, u1] = seq();
        pre::Vec3f point = tri.area_sample(u0, u1).point;
        if (bound.contains(point)) {
          if (mode == Mode::Normal)
            add(point, normal);
          else
            add(point, area);
        }
      }
    }
  }
  void write() {
    std::ofstream stream(
        outputFilename, std::ios_base::out | std::ios_base::binary);
    if (!stream.is_open()) {
      std::cerr << "Failed to open " + pre::show(outputFilename) << "!\n";
      std::exit(EXIT_FAILURE);
    }
    stream.write(reinterpret_cast<const char*>(count.data()), 8 * 3);
    uint64_t dim = 1;
    if (mode == Mode::Normal) dim = 3;
    stream.write(reinterpret_cast<const char*>(&dim), 8);
    stream.write(reinterpret_cast<const char*>(bound[0].data()), 4 * 3);
    stream.write(reinterpret_cast<const char*>(bound[1].data()), 4 * 3);
    stream.write(
        reinterpret_cast<const char*>(voxelData.data()),
        sizeof(voxelData[0]) * voxelData.size());
  }

  Mode mode = Mode::Area;
  pre::Bound3f bound;
  float voxelSizeZ = 0;
  pre::Vec3<float> voxelSize;
  pre::Vec3<float> invVoxelSize;
  pre::Vec3<uint64_t> count = {128, 128, 128};
  std::vector<float> voxelData;
  std::string outputFilename = "Voxels.bin";
  pre::monte_carlo::LowDiscrepancySequence2f seq;
  std::string includedIds;
  std::string excludedIds;
};

void Voxelize(const std::string& filename, VoxelRegion& region);

int main(int argc, char** argv) {
  programName = argv[0];
  std::string filename;
  std::string mode;
  VoxelRegion region;
  {
    pre::terminal::OptionParser options("[OPTIONS] FILENAME");
    options.on_option(
        "-s/--subdivs", &region.count[0], &region.count[1], &region.count[2]) =
        "Specify the number of subdivisions in X, Y, and Z.";
    options.on_option(
        "--lower",  //
        &region.bound[0][0], &region.bound[0][1], &region.bound[0][2]) =
        "Specify the minimum X, Y, and Z coordinates of the region to "
        "voxelize.";
    options.on_option(
        "--upper",  //
        &region.bound[1][0], &region.bound[1][1], &region.bound[1][2]) =
        "Specify the maximum X, Y, and Z coordinates of the region to "
        "voxelize.";
    options.on_option("-m/--mode", &mode) =
        "Specify the collection mode. Valid values are 'area', 'cosine-area', "
        "and 'normal'.";
    options.on_option("-z/--fit-z", &region.voxelSizeZ) =
        "Specify the desired Z voxel size. This lets the program fit the Z "
        "coordinate range to the scene geometry automatically.";
    options.on_option("-o/--output", &region.outputFilename) =
        "Specify the output filename, \"Voxels.bin\" by default.";
    options.on_option("--include-ids", &region.includedIds) =
        "Specify an ECMA regular expression matching the material IDs to "
        "include.";
    options.on_option("--exclude-ids", &region.excludedIds) =
        "Specify an ECMA regular expression matching the material IDs to "
        "exclude.";
    options.on_option("-y/--yes", &alwaysYes) =
        "Answer yes to all confirmation questions. Only use this if you know "
        "what you are doing!";
    options.on_option("--version", [] {
      std::cout << "Version 1.4" << std::endl;
      std::exit(EXIT_SUCCESS);
    }) = "Show the current version number.";
    options.on_positional([&filename](std::string_view argv) {
      if (not filename.empty())
        throw std::runtime_error("Expected only 1 filename!");
      filename = argv;
    });
    options.parse_or_exit(argc, argv);
  }
  if (filename.empty()) {
    std::cerr << "Expected a filename!\n";
    QuitWithHint();
  }
  if (!mode.empty()) {
    if (mode == "area")
      region.mode = VoxelRegion::Mode::Area;
    else if (mode == "cosine-area")
      region.mode = VoxelRegion::Mode::CosineArea;
    else if (mode == "normal")
      region.mode = VoxelRegion::Mode::Normal;
    else {
      std::cerr << "The mode is invalid!\n";
      std::cerr << "This must be 'area', 'cosine-area', or 'normal'.\n\n";
      QuitWithHint();
    }
  }
  if (!region.bound) {
    std::cerr << "The region is invalid!\nYou must specify the lower "
                 "and upper bounds of the region with '--lower' "
                 "and '--upper'.\n\n";
    QuitWithHint();
  }
  try {
    Voxelize(filename, region);
    return EXIT_SUCCESS;
  } catch (const std::exception& exception) {
    std::cerr << "Unhandled exception!\n";
    std::cerr << "exception.what(): " << exception.what() << '\n';
  } catch (const H5::Exception& exception) {
    std::cerr << "Unhandled exception from HDF5!";
    std::cerr << "\n  " << exception.getFuncName();
    std::cerr << "\n  " << exception.getDetailMsg();
  }
  return EXIT_FAILURE;
}

struct HDFData {
  struct MasterTableRecord {
    unsigned type = 0;
    unsigned typeIndex = 0;
    unsigned parentIndex = 0;
  };

  using MasterTableRecords = std::vector<MasterTableRecord>;

  static MasterTableRecords ReadMasterTable(H5::H5File& file) {
    MasterTableRecords records;
    auto dataSet = file.openDataSet("MasterTable");
    records.resize(dataSet.getSpace().getSimpleExtentNpoints());
    dataSet.read(records.data(), dataSet.getDataType());
    return records;
  }

  using MaterialNames = std::vector<std::string>;

  static MaterialNames ReadMaterialNames(H5::H5File& file) {
    MaterialNames names;
    auto dataSet = file.openGroup("Materials").openDataSet("Names");
    size_t dataSize = dataSet.getSpace().getSimpleExtentNpoints();
    std::vector<char*> cstrs(dataSize);
    static const H5::StrType cstrType(0, H5T_VARIABLE);
    dataSet.read(cstrs.data(), cstrType);
    names.resize(cstrs.size());
    for (size_t i = 0; i < cstrs.size(); i++) {
      std::string_view strv = cstrs[i];
      while (not strv.empty() and not strv.starts_with("(id = "))
        strv.remove_prefix(1);
      if (strv.size() < 7) continue;
      strv.remove_prefix(6);  // (id =
      strv.remove_suffix(1);  // )
      names[i] = std::string(strv);
    }
    if (!cstrs.empty())
      if (H5Dvlen_reclaim(
              cstrType.getId(), dataSet.getSpace().getId(), H5P_DEFAULT,
              cstrs.data()) < 0)
        throw std::runtime_error("Bad reclaim!");
    return names;
  }

  struct MaterialRemap {
    std::vector<std::pair<uint32_t, uint32_t>> ranges;
    std::vector<uint32_t> data;
  };

  static MaterialRemap ReadMaterialRemap(H5::H5File& file) {
    MaterialRemap remap;
    auto remapGroup = file.openGroup("Instances")
                          .openGroup("NoMotion")
                          .openGroup("MaterialRemapping");
    auto dataDataSet = remapGroup.openDataSet("Data");
    remap.data.resize(dataDataSet.getSpace().getSimpleExtentNpoints());
    dataDataSet.read(remap.data.data(), dataDataSet.getDataType());
    auto rangesDataSet = remapGroup.openDataSet("Ranges");
    remap.ranges.resize(rangesDataSet.getSpace().getSimpleExtentNpoints());
    rangesDataSet.read(remap.ranges.data(), rangesDataSet.getDataType());
    return remap;
  }

  struct FacetObject {
    pre::Bound3f bound;
    std::vector<pre::Vec3<float>> verts;
    std::vector<pre::Vec3<unsigned>> faces;
    std::vector<unsigned> materials;
    void read(H5::Group& group) {
      group.openDataSet("BoundingBox")
          .read(bound.data(), H5::PredType::NATIVE_FLOAT);
      auto vertsDataSet = group.openDataSet("Vertices");
      auto facesDataSet = group.openDataSet("Facets");
      size_t vertsDataSize = vertsDataSet.getSpace().getSimpleExtentNpoints();
      size_t facesDataSize = facesDataSet.getSpace().getSimpleExtentNpoints();
      THROW_IF(vertsDataSize % 3, "Expected all triangles!");
      THROW_IF(facesDataSize % 3, "Expected all triangles!");
      verts.resize(vertsDataSize / 3);
      faces.resize(facesDataSize / 3);
      vertsDataSet.read(verts.data(), H5::PredType::NATIVE_FLOAT);
      facesDataSet.read(faces.data(), H5::PredType::NATIVE_UINT);
      materials.resize(faces.size());
      group.openDataSet("FacetMaterials")
          .read(materials.data(), H5::PredType::NATIVE_UINT);
    }
    void Voxelize(
        const pre::Mat4f& transform,
        VoxelRegion& region,
        const uint32_t* remap,
        const std::vector<bool>& materialsOk) {
      pre::Mat3f linear = transform;
      pre::Vec3f affine = transform.col(3).head<3>();
      auto materialIter = materials.begin();
      for (const auto& face : faces) {
        if (materialsOk.empty() or materialsOk[remap[*materialIter++]]) {
          pre::Triangle3f tri = {
              pre::dot(linear, verts[face[0]]) + affine,
              pre::dot(linear, verts[face[1]]) + affine,
              pre::dot(linear, verts[face[2]]) + affine};
          region.add(tri);
        }
      }
    }
  };

  using FacetObjects = std::vector<FacetObject>;

  static FacetObjects ReadFacetObjects(H5::H5File& file) {
    auto group = file.openGroup("Objects").openGroup("FacetObjects");
    FacetObjects result;
    result.resize(group.getNumObjs());
    for (size_t index = 0; index < result.size(); index++) {
      std::string name = group.getObjnameByIdx(index);
      auto group2 = group.openGroup(name);
      result[std::stoi(name)].read(group2);
    }
    return result;
  }

  using Affine = pre::Array<float, 3, 4>;
  using Affines = std::vector<Affine>;
  static_assert(sizeof(Affine) == 12 * sizeof(float));

  static Affines ReadInstanceTransforms(H5::H5File& file) {
    Affines affines;
    auto dataSet = file.openGroup("Instances")
                       .openGroup("NoMotion")
                       .openDataSet("Transforms");
    affines.resize(dataSet.getSpace().getSimpleExtentNpoints());
    dataSet.read(affines.data(), dataSet.getDataType());
    return affines;
  }

  MasterTableRecords master;
  MaterialNames materialNames;
  MaterialRemap materialRemap;
  FacetObjects objects;
  Affines transforms;

  HDFData(const std::string& filename) {
    H5::H5File file(filename, H5F_ACC_RDONLY);
    master = ReadMasterTable(file);
    materialNames = ReadMaterialNames(file);
    materialRemap = ReadMaterialRemap(file);
    objects = ReadFacetObjects(file);
    transforms = ReadInstanceTransforms(file);
  }
};

struct Record {
  Record* parent = nullptr;
  HDFData::FacetObject* object = nullptr;
  pre::Mat4f transform = pre::Mat4f::identity();
  bool isLeaf = true;
  unsigned typeIndex = 0;

  HDFData::FacetObject* RootObject() {
    Record* record = this;
    while (record->parent) record = record->parent;
    return record->object;
  }

  pre::Mat4f FullTransform() const {
    return parent ? pre::dot(parent->FullTransform(), transform) : transform;
  }

  pre::Bound3f ObjectBound() const {
    THROW_IF(not object, "Object should be non-NULL by now!");
    pre::Bound3f bound = object->bound;
    pre::Bound3f result;
    for (auto [i, j, k] : pre::ArrayRange(2, 2, 2))
      result |= pre::homogeneous_dot(
          transform, pre::Vec3f{bound[i][0], bound[j][1], bound[k][2]});
    return result;
  }

  void TrueZBound(const VoxelRegion& region, float& zmin, float& zmax) const {
    THROW_IF(not object, "Object should be non-NULL by now!");
    pre::Mat3f linear = transform;
    pre::Vec3f affine = transform.col(3).head<3>();
    for (pre::Vec3f v : object->verts) {
      v = pre::dot(linear, v) + affine;
      if (pre::Bound2f(region.bound).contains(pre::Vec2f(v))) {
        zmin = std::min(zmin, v[2]);
        zmax = std::max(zmax, v[2]);
      }
    }
  }
};

using Records = std::vector<Record>;

void Voxelize(const std::string& filename, VoxelRegion& region) {
  HDFData hdfData(filename);
  std::vector<bool> materialsOk;
  if (not region.includedIds.empty() or  //
      not region.excludedIds.empty()) {
    materialsOk.resize(hdfData.materialNames.size(), true);
    bool anyExcluded = false;
    if (not region.includedIds.empty()) try {  // Exclude non-matches
        std::regex regex(region.includedIds);
        for (size_t i = 0; i < hdfData.materialNames.size(); i++)
          if (not std::regex_match(hdfData.materialNames[i], regex))
            materialsOk[i] = false, anyExcluded = true;
      } catch (const std::regex_error& error) {
        std::cerr << "Exception while handling regex "
                  << pre::show(region.includedIds) << '\n';
        std::cerr << "exception.what(): " << error.what() << '\n';
        std::exit(EXIT_FAILURE);
      }
    if (not region.excludedIds.empty()) try {  // Exclude matches
        std::regex regex(region.excludedIds);
        for (size_t i = 0; i < hdfData.materialNames.size(); i++)
          if (std::regex_match(hdfData.materialNames[i], regex))
            materialsOk[i] = false, anyExcluded = true;
      } catch (const std::regex_error& error) {
        std::cerr << "Exception while handling regex "
                  << pre::show(region.excludedIds) << '\n';
        std::cerr << "exception.what(): " << error.what() << '\n';
        std::exit(EXIT_FAILURE);
      }

    int totalOk = 0;
    for (bool ok : materialsOk)
      if (ok) totalOk++;
    std::cerr << "Excluded " << materialsOk.size() - totalOk << " of "
              << materialsOk.size() << " materials from the voxelization."
              << std::endl;
    if (anyExcluded) {
      std::cerr << "Included materials after filtering:" << std::endl;
      for (size_t i = 0; i < hdfData.materialNames.size(); i++)
        if (materialsOk[i])
          std::cerr << "  " << hdfData.materialNames[i] << std::endl;
    }
  }
  Records records(hdfData.master.size());
  for (size_t index = 0; index < records.size(); index++) {
    auto& record = records[index];
    unsigned type = hdfData.master[index].type;
    unsigned typeIndex = hdfData.master[index].typeIndex;
    unsigned parentIndex = hdfData.master[index].parentIndex;
    if (parentIndex < records.size()) {
      record.parent = &records[parentIndex];
      record.parent->isLeaf = false;
    }
    record.typeIndex = typeIndex;
    if (type != 1) record.isLeaf = false;
    if (type == 1) {  // Static instance?
      record.transform = hdfData.transforms[typeIndex];
      record.transform[3][3] = 1;
    } else if (type == 16) {  // Facet object?
      record.object = &hdfData.objects[typeIndex];
    }
  }
  float zmin = +INFINITY;
  float zmax = -INFINITY;
  std::vector<Record*> recordsToVoxelize;
  for (auto& record : records)
    if (record.isLeaf) {
      auto* object = record.RootObject();
      if (not object) continue;
      record.object = object;
      record.transform = record.FullTransform();
      pre::Bound3f bound = record.ObjectBound();
      if (region.voxelSizeZ > 0) {
        if (pre::Bound2f(region.bound).overlaps(bound)) {
          record.TrueZBound(region, zmin, zmax);
          recordsToVoxelize.push_back(&record);
        }
      } else if (region.bound.overlaps(bound)) {
        recordsToVoxelize.push_back(&record);
      }
    }
  if (recordsToVoxelize.empty()) {
    std::cerr << "There is no geometry to voxelize in the specified region!\n"
                 "HINT: If you use the option '--fit-z', the program will fit "
                 "the Z coordinates of the region to the scene geometry.\n\n";
    QuitWithHint();
  }
  if (region.voxelSizeZ > 0) {
    region.bound[0][2] = zmin;
    region.bound[1][2] = zmax;
    region.count[2] = pre::max(1.f, (zmax - zmin) / region.voxelSizeZ);
    std::cerr << "Z min = " << zmin << " meters\n";
    std::cerr << "Z max = " << zmax << " meters\n";
    std::cerr << "Setting Z subdivisions to " << region.count[2];
    std::cerr << "...\n";
  }
  double GiB =
      static_cast<long double>(region.count.prod() * 4ULL) / 1073741824.0L;
  if (GiB > 1 and not alwaysYes) {
    std::stringstream message;
    message << std::fixed;
    message << std::setprecision(2);
    message << "WARNING: The resulting dataset is going to be " << GiB;
    message << " GiB on disk. Are you sure you want to continue?";
    if (!pre::terminal::ask_yes_no(message.str(), false)) {
      std::cerr << "Terminating...\n";
      std::exit(EXIT_SUCCESS);
    }
  }
  region.invVoxelSize = region.count / region.bound.extent();
  region.voxelSize = region.bound.extent() / region.count;
  if (region.mode == VoxelRegion::Mode::Normal)
    region.voxelData.resize(region.count.prod() * 3);
  else
    region.voxelData.resize(region.count.prod());
  int i = 1;
  int n = recordsToVoxelize.size();
  for (auto* record : recordsToVoxelize) {
    std::cerr << "\r";
    std::cerr << "Voxelizing [" << i++ << "/" << n << "]";
    record->object->Voxelize(
        record->transform, region,
        hdfData.materialRemap.data.data() +
            hdfData.materialRemap.ranges[record->typeIndex].first,
        materialsOk);
  }
  std::cerr << std::endl;
  if (region.mode == VoxelRegion::Mode::Normal) {
    for (uint64_t i = 0; i < region.count.prod(); i++) {
      float denom = 1.0f / pre::length(pre::Vec3f(&region.voxelData[i * 3]));
      if (std::isfinite(denom)) {
        region.voxelData[i * 3 + 0] *= denom;
        region.voxelData[i * 3 + 1] *= denom;
        region.voxelData[i * 3 + 2] *= denom;
      } else {
        region.voxelData[i * 3 + 0] = 0;
        region.voxelData[i * 3 + 1] = 0;
        region.voxelData[i * 3 + 2] = 0;
      }
    }
  }
  region.write();
}
