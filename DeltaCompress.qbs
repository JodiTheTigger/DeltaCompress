import qbs

CppApplication {
    type: "application" // To suppress bundle generation on Mac
    consoleApplication: true

    Depends { name: "cpp" }

    // http://stackoverflow.com/questions/25936548/how-to-use-qbs-to-build-a-qt-application-with-qt-creator
    //cpp.cxxLanguageVersion: "c++11"
    // Clang is annoying and warns missing braces -Wall.
    // Have a static array used as a default parameter and nowhere else
    // again, clang is annoying and complains.
    cpp.cxxFlags:
    [
        "-std=c++14",
        "-Wextra",
        "-Wno-missing-braces",
        "-Wno-unneeded-internal-declaration"
    ]
    cpp.treatWarningsAsErrors: true
    cpp.warningLevel: "all"

    files: [
        "DeltaCompress.cpp",
        "Range_coding.hpp",
    ]

    Group {     // Properties for the produced executable
        fileTagsFilter: product.type
        qbs.install: true
    }
}

