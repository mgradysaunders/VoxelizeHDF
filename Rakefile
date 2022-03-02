require 'rake'

CXX = "g++-11"
CMAKE = "cmake"
CMAKE_PREFIX_PATH = ["/home/grady/git/mgradysaunders/Precept/build/install",
                     "/home/grady/git/HDFGroup/hdf5/build/install"]
CMAKE_FLAGS = "-S. -Bbuild -GNinja -DCMAKE_CXX_COMPILER=\"#{CXX}\" -DCMAKE_PREFIX_PATH=\"#{CMAKE_PREFIX_PATH.join ';'}\" -DCMAKE_INSTALL_PREFIX=build/install"

desc "Default task."
task :default do
  # nothing
end

desc "Build executable."
task :build do
  sh "#{CMAKE} #{CMAKE_FLAGS}" unless File.directory? 'build'
  sh "#{CMAKE} --build build"
  sh "rm -f compile_commands.json"
  sh "ln -s build/compile_commands.json ."
end

task :build_image do
  sh "rm -rf VoxelizeHDF-x86_64.AppDir"
  sh "linuxdeploy --appdir=VoxelizeHDF-x86_64.AppDir -e build/VoxelizeHDF -d VoxelizeHDF.desktop -i icon.svg"
  sh "appimagetool VoxelizeHDF-x86_64.AppDir/ VoxelizeHDF-x86_64.AppImage"
  sh "rm -rf VoxelizeHDF-x86_64.AppDir"
end
