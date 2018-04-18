//
// Developed by:  <Abd Elrahman Bachet> (abdelrahman.bachet@gmail.com)
//                http://www.biodataanalysis.de/
//
// With contributions by:
//
//
// Copyright (c) 2014-2018, BioDataAnalysis GmbH
// All Rights Reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are not permitted. All information contained herein
// is, and remains the property of BioDataAnalysis GmbH.
// Dissemination of this information or reproduction of this material
// is strictly forbidden unless prior written permission is obtained
// from BioDataAnalysis GmbH.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//

#include <vigra/imageinfo.hxx>
#include <vigra/impex.hxx>
#include <vigra/multi_array.hxx>
#include <vigra/stdimage.hxx>

#include <iostream>
#include <string>

using namespace vigra;

int main(int argc, char** argv) {
    std::cout << "Task 01 started" << std::endl;

    // You can add your code here

    /* The task is to crop an image and save the output image with a new Name
    Helpful links: https://ukoethe.github.io/vigra/doc-release/vigra/ImageInputOutputTutorial.html */

    //load the image
    const char* inputImage;
    inputImage = "../images/bDZ17-1I_wE02_s7_z1_t1_cGFP_u001.tif";
    ImageImportInfo imageInfo(inputImage, 0);

    //copy image pixels into multiarray data structure
    MultiArray<2, UInt8> imageOriginal(imageInfo.shape());
    importImage(imageInfo, imageOriginal);

    //create a new array for output image of size 500x350
    MultiArray<2, UInt8> imageCropped(Shape2(500, 350));

    //crop the image starting from index 400,500
    for (int i = 0; i < imageCropped.shape(1); ++i) {
        for (int j = 0; j < imageCropped.shape(0); ++j) {
            imageCropped(j, i) = imageOriginal(400 + j, 500 + i);
        }
    }

    //save the output image
    exportImage(imageCropped, "../images/croppedImage.tif");

    std::cout << "Task 01 finished successfully" << std::endl;
    return 0;
}
