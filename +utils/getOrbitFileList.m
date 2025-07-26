function fileList = getOrbitFileList()
    files = dir('filtered PlanarOrbitData/*.csv');
    fileList = {files.name};
end