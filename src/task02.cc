//
// Developed by:  <Your Name> (your@email)
//                http://www.biodataanalysis.de/
//
// With contributions by:
//
//
// Copyright (c) 2024-2028, BioDataAnalysis GmbH
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

#include <vigra/multi_array.hxx>
#include <vigra/impex.hxx>
#include <vigra/stdconvolution.hxx>
#include <vigra/convolution.hxx>

#include <iostream>
#include <string>
#include <cmath>
#include <vector>

using namespace vigra;

typedef std::vector<std::vector <double>> Matrix;

Matrix gaussianFilter(int kernelSize, double sigma) {
    Matrix kernel(kernelSize, std::vector<double>(kernelSize));
    double sum = 0.0;

    for (int i= -((kernelSize-1)/2) ; i<(kernelSize+1)/2 ; i++) {
        for (int j=-((kernelSize-1)/2) ; j<(kernelSize+1)/2 ; j++) {
            kernel[i+((kernelSize-1)/2)][j+((kernelSize-1)/2)] = exp(-(i*i+j*j)/(2*sigma*sigma))/(2*M_PI*sigma*sigma);
            sum += kernel[i+((kernelSize-1)/2)][j+((kernelSize-1)/2)];
        }
    }
    //to check:
    //std::cout<<"sum = "<<sum<<std::endl;

    for (int i= -((kernelSize-1)/2) ; i<(kernelSize+1)/2 ; i++) {
        for (int j=-((kernelSize-1)/2) ; j<(kernelSize+1)/2 ; j++) {
            kernel[i+((kernelSize-1)/2)][j+((kernelSize-1)/2)] /= sum;
        }
    }

    return kernel;
}

int main(int argc, char** argv)
{
    std::cout << "Task 02 started" << std::endl;

    // You can add your code here

    /* The Task is to smooth an image using a Gaussian filter
    Helpful links: https://ukoethe.github.io/vigra/doc-release/vigra/ImageProcessingTutorial.html#SmoothingTutorial
                   https://homepages.inf.ed.ac.uk/rbf/HIPR2/gsmooth.htm */

    //load the image
    const char* inputImage;
    inputImage = "../images/bDZ17-1I_wE02_s7_z1_t1_cDAPI_u001.tif";
    ImageImportInfo imageInfo(inputImage,0);

    //copy image pixels into multiarray data structure
    MultiArray<2, UInt8> imageOriginal(imageInfo.shape());
    importImage(imageInfo, imageOriginal);

    //set the gaussian filter parameters
    int kernelSize = 11;
    double sigma = 3.0;
    int right = (kernelSize-1)/2;
    int left = -right;

    //apply the filter
    Kernel2D<double> filter;
    filter.initExplicitly(Shape2(left,left), Shape2(right,right));
    filter.initGaussian(sigma);
    MultiArray<2, UInt8> imageSmoothed(imageInfo.shape());

    //another way to generate a gaussianFilter manually

    /* Matrix gFilter = gaussianFilter(kernelSize,sigma);

    for (int i=left;i<=right;++i){
          for (int j=left;i<=right;++j){
              filter(i,j) = gFilter[i-left][j-left];
              std::cout<<filter(i,j);
      }
      std::cout<<std::endl;
    }*/

    //convolve the image
    convolveImage(imageOriginal, imageSmoothed, filter);

    //save the output image
    exportImage(imageSmoothed, "../images/smoothedImage.tif");


    std::cout << "Task 02 finished successfully" << std::endl;
    return 0;
}
