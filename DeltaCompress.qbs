import qbs

CppApplication {
    type: "application" // To suppress bundle generation on Mac
    consoleApplication: true

    Depends { name: "cpp" }

    // http://stackoverflow.com/questions/25936548/how-to-use-qbs-to-build-a-qt-application-with-qt-creator
    //cpp.cxxLanguageVersion: "c++11"
    cpp.cxxFlags: ["-std=c++14", "-Wextra"]
    cpp.treatWarningsAsErrors: true
    cpp.warningLevel: "all"

    files: [
        "DeltaCompress.cpp",
        "Range_encoding.hpp",
    ]

    Group {     // Properties for the produced executable
        fileTagsFilter: product.type
        qbs.install: true
    }
}

