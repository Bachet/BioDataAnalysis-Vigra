//
// Developed by:  <Abd Elrahman Bachet> (abdelrahman.bachet@gmail.com)
//                http://www.biodataanalysis.de/
//
// With contributions by:
//
//
// Copyright (c) 2034-2038, BioDataAnalysis GmbH
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

#include <vigra/impex.hxx>
#include <vigra/multi_array.hxx>

#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>

using namespace vigra;

//macros
#define bitDepth 8

void histogramValues(long int bitsPerPixel, int* histogram) {
    std::ofstream fout("../images/histogram.txt");
    for (long int i = 0; i < bitsPerPixel; ++i) {
        fout << i << " : " << histogram[i] << std::endl;
    }
}

long int otsu(long int bitsPerPixel, int* histogram, int64_t N) {

    int64_t q1 = 0, q2 = 0, m1 = 0, m2 = 0, sum = 0, sum2 = 0;
    long int threshold;
    double betweenVar, maxVar = 0;

    for (long int i = 0; i < bitsPerPixel; ++i) {
        sum += i * histogram[i];
    }

    for (long int i = 0; i < bitsPerPixel; ++i) {
        q1 += histogram[i];
        q2 = N - q1;

        sum2 += i * histogram[i];
        if (q1 != 0 && q2 != 0) {
            m1 = sum2 / q1;
            m2 = (sum - sum2) / q2;
        } else {
            continue;
        }

        betweenVar = q1 * q2 * pow((m1 - m2), 2);

        //maximizing variance between the two regions
        if (betweenVar > maxVar) {
            threshold = i;
            maxVar = betweenVar;
        }
    }

    return threshold;
}

int main(int argc, char** argv) {
    std::cout << "Task 03 started" << std::endl;

    // You can add your code here

    /* The Task is to apply image segmentation using OTSU's method to an image
    Helpful links: http://www.math.tau.ac.il/~turkel/notes/otsu.pdf
                   http://www.ipol.im/pub/art/2016/158/ */

    //load the image
    const char* inputImage;
    inputImage = "../images/bDZ17-1I_wE02_s7_z1_t1_cDAPI_u001.tif";
    ImageImportInfo imageInfo(inputImage, 0);

    //copy image pixels into multiarray data structure
    MultiArray<2, UInt8> imageOriginal(imageInfo.shape());
    importImage(imageInfo, imageOriginal);

    long int bitsPerPixel = pow(2, bitDepth);
    int* histogram = new int[bitsPerPixel];
    for (long int i = 0; i < bitsPerPixel; ++i) {
        histogram[i] = 0;
    }

    clock_t begin = clock();
    for (int i = 0; i < imageInfo.height(); ++i) {
        for (int j = 0; j < imageInfo.width(); ++j) {
            histogram[imageOriginal(j, i)]++;
        }
    }

    //for insights: writing histogram output into file
    //histogramValues(bitsPerPixel, histogram);

    //applying otsu algorithm
    int64_t N = imageInfo.height() * imageInfo.width();

    long int threshold = otsu(bitsPerPixel, histogram, N);
    clock_t end = clock();

    std::cout << "threshold is: " << threshold << std::endl;

    //calculate time taken to implement the otsu algorithm (o(N))
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "Time taken for implementation: " << elapsed_secs << "s" << std::endl;

    MultiArray<2, UInt8> outputImage(imageInfo.shape());

    for (int i = 0; i < imageInfo.height(); ++i) {
        for (int j = 0; j < imageInfo.width(); ++j) {
            if (imageOriginal(j, i) > threshold) {
                outputImage(j, i) = bitsPerPixel - 1;
            } else {
                outputImage(j, i) = 0;
            }
        }
    }

    //save the output image
    exportImage(outputImage, "../images/segmentedImage.tif");

    std::cout << "Task 03 finished successfully" << std::endl;
    return 0;
}
