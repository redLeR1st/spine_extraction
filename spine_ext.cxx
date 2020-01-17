/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include <iostream>
#include <map>
#include <set>
#include <algorithm>
#include <functional>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"
#include "itkPasteImageFilter.h"

#include "itkBinaryThresholdImageFilter.h"

 // Software Guide : BeginCodeSnippet
#include "itkHoughTransform2DCirclesImageFilter.h"
 // Software Guide : EndCodeSnippet
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkThresholdImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"
#include <list>
#include "itkCastImageFilter.h"
#include "itkMath.h"

#include "itkFlipImageFilter.h"

#include "itkTileImageFilter.h"

#include <itkSpatialObject.h>

#include <itkCropImageFilter.h>

#include "itkScaleTransform.h"

#include "itkResampleImageFilter.h"
#include "itkVersor.h"
#include "itkChangeInformationImageFilter.h"
#include "itkImageDuplicator.h"

#include "itkOtsuMultipleThresholdsImageFilter.h"

#include "itkEllipseSpatialObject.h"

#include "itkBinaryMorphologicalClosingImageFilter.h"
#include "itkBinaryBallStructuringElement.h"

#include "itkScalarImageToHistogramGenerator.h"
#include <itkOtsuMultipleThresholdsCalculator.h>

using PixelType = short;

using ImageType = itk::Image< PixelType, 3 >;
using ImageType2D = itk::Image< PixelType, 2 >;

using ReaderType = itk::ImageFileReader< ImageType >;
using WriterType = itk::ImageFileWriter< ImageType >;

typedef itk::SpatialObject<2>             SpatialObjectType;
typedef SpatialObjectType::TransformType  TransformType;
typedef   float           AccumulatorPixelType;

typedef itk::EllipseSpatialObject<2> ElipseType;

typedef itk::HoughTransform2DCirclesImageFilter<PixelType,
    AccumulatorPixelType> HoughTransformFilterType;
typedef HoughTransformFilterType::CirclesListType CirclesListType;
//typedef HoughTransformFilterType::CircleType::Type CirclesListType;

typedef itk::OtsuMultipleThresholdsImageFilter <ImageType2D, ImageType2D> OtsuFilterType2D;
typedef itk::OtsuMultipleThresholdsImageFilter <ImageType, ImageType> OtsuFilterType;

using ScalarImageToHistogramGeneratorType = itk::Statistics::ScalarImageToHistogramGenerator<ImageType>;
using HistogramType = ScalarImageToHistogramGeneratorType::HistogramType;

using CalculatorType = itk::OtsuMultipleThresholdsCalculator<HistogramType>;


bool circle_itersect(double x1, double y1, double r1, double x2, double y2, double r2) {
    double distSq = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
    double radSumSq_m = (r1 - r2) * (r1 - r2);
    double radSumSq_p = (r1 + r2) * (r1 + r2);
    
    if (radSumSq_m <= distSq && distSq <= radSumSq_p) { // intersect
        return true;
    }
    else     // not intersect
        return false;
}

void corrigate_circles(std::vector<CirclesListType> &circles_list) {
    
    TransformType::OffsetType Object2ToObject1Offset;
    Object2ToObject1Offset[0] = 100;
    Object2ToObject1Offset[1] = 100;

    double av_x = 0;
    double av_y = 0;
    double av_rad = 0;

    int sum = 0;

    auto it = circles_list.begin();
    while (it != circles_list.end())
    {

        CirclesListType::const_iterator itCircles = (*it).begin();
        while (itCircles != (*it).end())
        {
           
            
            av_x += (*itCircles)->GetObjectToParentTransform()->GetOffset()[0];
            av_y += (*itCircles)->GetObjectToParentTransform()->GetOffset()[1];
            av_rad += (*itCircles)->GetRadius()[0];

            sum++;
            itCircles++;
        }
        it++;
    }
    av_x /= sum;
    av_y /= sum;
    av_rad /= sum;

    std::cerr << "AVG_X " << av_x << std::endl;
    std::cerr << "AVG_Y " << av_y << std::endl;
    std::cerr << "AVG_RAD " << av_rad << std::endl;

    int filter_size = 40;

    it = circles_list.begin();
    while (it != circles_list.end())
    {
        
        CirclesListType::const_iterator itCircles = (*it).begin();
        
        //std::cout << "x1 " << (*itCircles)->GetObjectToParentTransform()->GetOffset()[0] << std::endl;
        //std::cout << "y1 " << (*itCircles)->GetObjectToParentTransform()->GetOffset()[1] << std::endl;
        //std::cout << "r1 " << (*itCircles)->GetRadius()[0]<< std::endl;
        //
        //std::cout << "x2 " << (*(*it2).begin())->GetObjectToParentTransform()->GetOffset()[0] << std::endl;
        //std::cout << "y2 " << (*(*it2).begin())->GetObjectToParentTransform()->GetOffset()[1] << std::endl;
        //std::cout << "r2 " << (*(*it2).begin())->GetRadius()[0] << std::endl;
        //std::cout << "INTER " << circle_itersect(av_x, av_y, av_rad,26,172,10) << std::endl;
        //std::cout << "NOT INTER " << circle_itersect(3,4,0,14,18,0) << std::endl;
        
        
        //     if (!circle_itersect(av_x, av_y, av_rad, (*(*it2).begin())->GetObjectToParentTransform()->GetOffset()[0],
        //         (*(*it2).begin())->GetObjectToParentTransform()->GetOffset()[1], (*(*it2).begin())->GetRadius()[0])) {
        if (!circle_itersect((*itCircles)->GetObjectToParentTransform()->GetOffset()[0], (*itCircles)->GetObjectToParentTransform()->GetOffset()[1], (*itCircles)->GetRadius()[0],
            av_x, av_y, av_rad)) {

            std::cout << "CORRIGATE NEEDED " << std::endl;

            bool has_inter = false;
            int min_dist = 99999999;
            while (++itCircles != (*it).end()) {
               int x = (*itCircles)->GetObjectToParentTransform()->GetOffset()[0];
               int y = (*itCircles)->GetObjectToParentTransform()->GetOffset()[1];
                double dist = (x - av_x) * (x - av_x) + (y - av_y) * (y - av_y);
                
                if (min_dist > dist) {

                    CirclesListType::const_iterator first = (*it).begin();
                    
                    Object2ToObject1Offset[0] = (*itCircles)->GetObjectToParentTransform()->GetOffset()[0];
                    Object2ToObject1Offset[1] = (*itCircles)->GetObjectToParentTransform()->GetOffset()[1];
                    (*first)->GetObjectToParentTransform()->SetOffset(Object2ToObject1Offset);
                    (*first)->SetRadius((*itCircles)->GetRadius()[0]);
                    min_dist = dist;
                }

                //if (circle_itersect((*itCircles)->GetObjectToParentTransform()->GetOffset()[0], (*itCircles)->GetObjectToParentTransform()->GetOffset()[1], (*itCircles)->GetRadius()[0],
                //    av_x, av_y, av_rad)) {
                //    
                //    CirclesListType::const_iterator first = (*it).begin();
                //
                //    Object2ToObject1Offset[0] = (*itCircles)->GetObjectToParentTransform()->GetOffset()[0];
                //    Object2ToObject1Offset[1] = (*itCircles)->GetObjectToParentTransform()->GetOffset()[1];
                //    (*first)->GetObjectToParentTransform()->SetOffset(Object2ToObject1Offset);
                //    (*first)->SetRadius((*itCircles)->GetRadius()[0]);
                //    break;
                //    has_inter = true;
                //
                //}
            }
            
            if (!has_inter) {
                std::cout << "NO BETTER CIRCLE FOUND" << std::endl;
            }

            // MEDIAN METHOD!
            // ///////////////////////////////////////////
            // std::map<int, double> dis_from_avg;
            // int key = 0;
            // 
            // bool flag = false;
            // auto it2 = std::next(it, 0);
            // //for (int i = 0; i < filter_size; i++) {
            // while (key < filter_size) {
            //     if (++it2 == circles_list.end()) {
            //         flag = true;
            //         break;
            //     }
            //     else {
            //         int x = (*(*it2).begin())->GetObjectToParentTransform()->GetOffset()[0];
            //         int y = (*(*it2).begin())->GetObjectToParentTransform()->GetOffset()[1];
            //         int rad = (*(*it2).begin())->GetRadius()[0];
            // 
            //         if (circle_itersect(x, y, rad, av_x, av_y, av_rad)) {
            //             double dist = (x - av_x) * (x - av_x) + (y - av_y) * (y - av_y);
            //             dis_from_avg.insert(std::make_pair(key, dist));
            //             key++;
            //         }
            // 
            //     }
            // }
            // if (flag) {
            //     break;
            // }
            // typedef std::pair<int, double> pair;
            // // create a empty vector of pairs
            // std::vector<pair> vec;
            // 
            // // copy key-value pairs from the map to the vector
            // std::copy(dis_from_avg.begin(),
            //     dis_from_avg.end(),
            //     std::back_inserter<std::vector<pair>>(vec));
            // 
            // // sort the vector by increasing order of its pair's second value
            // // if second value are equal, order by the pair's first value
            // std::sort(vec.begin(), vec.end(),
            //     [](const pair& l, const pair& r) {
            //     if (l.second != r.second)
            //         return l.second < r.second;
            // 
            //     return l.first < r.first;
            // });
            // 
            // // print the vector
            // for (auto const &pair : vec) {
            //     std::cout << '{' << pair.first << "," << pair.second << '}' << '\n';
            // }
            // 
            // auto median_key_iterator = std::next(vec.begin(), filter_size / 3);
            // 
            // int median_circle_key = (*median_key_iterator).first;
            // 
            // auto median_circle = std::next(it, median_circle_key);
            // 
            // CirclesListType::const_iterator itCircles = (*it).begin();
            // 
            // Object2ToObject1Offset[0] = (*(*median_circle).begin())->GetObjectToParentTransform()->GetOffset()[0];
            // Object2ToObject1Offset[1] = (*(*median_circle).begin())->GetObjectToParentTransform()->GetOffset()[1];
            // (*itCircles)->GetObjectToParentTransform()->SetOffset(Object2ToObject1Offset);
            // (*itCircles)->SetRadius((*(*median_circle).begin())->GetRadius()[0]);
            // 
            // 
            // ///////////////////////////////////////////
            // 
            // 
            // 
            // //do {
            // //
            // //    Object2ToObject1Offset[0] = (*(*it2).begin())->GetObjectToParentTransform()->GetOffset()[0];
            // //    Object2ToObject1Offset[1] = (*(*it2).begin())->GetObjectToParentTransform()->GetOffset()[1];
            // //    (*itCircles)->GetObjectToParentTransform()->SetOffset(Object2ToObject1Offset);
            // //    (*itCircles)->SetRadius((*(*it2).begin())->GetRadius()[0]);
            // //
            // //    if (++it2 == circles_list.end()) {
            // //        break;
            // //    }
            // //}// while (circle_itersect(av_x, av_y, av_rad, (*(*it2).begin())->GetObjectToParentTransform()->GetOffset()[0],
            // // //   (*(*it2).begin())->GetObjectToParentTransform()->GetOffset()[1], (*(*it2).begin())->GetRadius()[0]));
            // //while (!circle_itersect((*itCircles)->GetObjectToParentTransform()->GetOffset()[0], (*itCircles)->GetObjectToParentTransform()->GetOffset()[1], (*itCircles)->GetRadius()[0], (*(*it2).begin())->GetObjectToParentTransform()->GetOffset()[0],
            // //    (*(*it2).begin())->GetObjectToParentTransform()->GetOffset()[1], (*(*it2).begin())->GetRadius()[0]));
            
        }
        //CirclesListType::const_iterator h = (*it).begin();
        // lets see where is the average
        //Object2ToObject1Offset[0] = av_x;
        //Object2ToObject1Offset[1] = av_y;
        //(*h)->GetObjectToParentTransform()->SetOffset(Object2ToObject1Offset);
        //(*h)->SetRadius(av_rad);
        it++;
        
    }

}

template<class MyType>
void draw_circle(MyType &input, std::vector<CirclesListType> circles_list) {

    ImageType::IndexType localIndex3D;

    const ImageType * inputImage = input->GetOutput();
    ImageType::RegionType inputRegion = inputImage->GetBufferedRegion();
    ImageType::SizeType size = inputRegion.GetSize();
    int tolerance = 100;
    int actual_slice_number = 0;
    
    std::cerr << "LIST SIZE " << circles_list.size() << std::endl;
    auto it = circles_list.begin();
    while (it != circles_list.end())
    {

        CirclesListType::const_iterator itCircles = (*it).begin();
        //while (itCircles != (*it).end())
        //{

            for (unsigned short x = 0; x < size[0]; x++) {
                for (unsigned short y = 0; y < size[1]; y++) {
                    int x_cube = (x - (*itCircles)->GetObjectToParentTransform()->GetOffset()[0])*(x - (*itCircles)->GetObjectToParentTransform()->GetOffset()[0]);
                    int y_cube = (y - (*itCircles)->GetObjectToParentTransform()->GetOffset()[1])*(y - (*itCircles)->GetObjectToParentTransform()->GetOffset()[1]);
                    
                    // keep only inside
                    if ((x_cube + y_cube) > ((*itCircles)->GetRadius()[0] + tolerance) * ((*itCircles)->GetRadius()[0]) + tolerance) {
                        localIndex3D[0] = x;
                        localIndex3D[1] = y;
                        localIndex3D[2] = actual_slice_number;
                        input->GetOutput()->SetPixel(localIndex3D, 0);
                        //input->GetOutput()->SetPixel(localIndex3D, -1);
                    }

                    // on the circle, then we draw
                    //tolerance = 0;
                    //int circle_line_tik = 65;
                    //if ((x_cube + y_cube) > ((*itCircles)->GetRadius()[0] + tolerance) * ((*itCircles)->GetRadius()[0]) + tolerance - circle_line_tik && (x_cube + y_cube) < ((*itCircles)->GetRadius()[0] + tolerance) * ((*itCircles)->GetRadius()[0]) + tolerance + circle_line_tik) {
                    //    localIndex3D[0] = x;
                    //    localIndex3D[1] = y;
                    //    localIndex3D[2] = actual_slice_number;
                    //    input->GetOutput()->SetPixel(localIndex3D, -10000);
                    //}

                }
            }
            //itCircles++;
        //}
        actual_slice_number++;
        it++;
    }
}


template<class MyType, class MyType2>
MyType do_hough_on_image(MyType input, MyType2 &original, bool keep_only_inside, int actual_slice_number, std::vector<CirclesListType> &circles_list) {
    //#########################################

#if 0
    if (argc < 6)
    {
        std::cerr << "Missing Parameters " << std::endl;
        std::cerr << "Usage: " << argv[0] << std::endl;
        std::cerr << " inputImage " << std::endl;      1
            std::cerr << " outputImage" << std::endl;      2
            std::cerr << " numberOfCircles " << std::endl; 3
            std::cerr << " radius Min " << std::endl;      4
            std::cerr << " radius Max " << std::endl;      5
            std::cerr << " sweep Angle (default = 0)" << std::endl; 6
            std::cerr << " SigmaGradient (default = 1) " << std::endl; 7
            std::cerr << " variance of the accumulator blurring (default = 5) " << std::endl; 8
            std::cerr << " radius of the disk to remove from the accumulator (default = 10) " << std::endl; 9
            return EXIT_FAILURE;
    }

#endif


    ImageType2D::IndexType localIndex;
    //ImageType::IndexType localIndex3D;
    typedef itk::Image< AccumulatorPixelType, 2 > AccumulatorImageType;

    ImageType2D::Pointer localImage = input->GetOutput();

    //  Software Guide : BeginLatex
    //
    //  We create the HoughTransform2DCirclesImageFilter based on the pixel
    //  type of the input image (the resulting image from the
    //  ThresholdImageFilter).
    //
    //  Software Guide : EndLatex
    // Software Guide : BeginCodeSnippet
    std::cout << "Computing Hough Map on slice: "<< actual_slice_number << std::endl;
    
    HoughTransformFilterType::Pointer houghFilter
        = HoughTransformFilterType::New();
    // Software Guide : EndCodeSnippet
    //  Software Guide : BeginLatex
    //
    //  We set the input of the filter to be the output of the
    //  ImageFileReader. We set also the number of circles we are looking for.
    //  Basically, the filter computes the Hough map, blurs it using a certain
    //  variance and finds maxima in the Hough map. After a maximum is found,
    //  the local neighborhood, a circle, is removed from the Hough map.
    //  SetDiscRadiusRatio() defines the radius of this disc proportional to
    //  the radius of the disc found.  The Hough map is computed by looking at
    //  the points above a certain threshold in the input image. Then, for each
    //  point, a Gaussian derivative function is computed to find the direction
    //  of the normal at that point. The standard deviation of the derivative
    //  function can be adjusted by SetSigmaGradient(). The accumulator is
    //  filled by drawing a line along the normal and the length of this line
    //  is defined by the minimum radius (SetMinimumRadius()) and the maximum
    //  radius (SetMaximumRadius()).  Moreover, a sweep angle can be defined by
    //  SetSweepAngle() (default 0.0) to increase the accuracy of detection.
    //
    //  The output of the filter is the accumulator.
    //
    //  Software Guide : EndLatex
    // Software Guide : BeginCodeSnippet

    houghFilter->SetInput(input->GetOutput());
    //houghFilter->SetNumberOfCircles(atoi(argv[3]));
    houghFilter->SetNumberOfCircles(10);

   houghFilter->SetMinimumRadius(10);
   houghFilter->SetMaximumRadius(30);
   houghFilter->SetSweepAngle(1);
   houghFilter->SetSigmaGradient(3);
   houghFilter->SetVariance(2);

    //houghFilter->SetMinimumRadius(18);
    ////houghFilter->SetMinimumRadius(10);
    //houghFilter->SetMaximumRadius(25);
    //houghFilter->SetSweepAngle(0);
    //houghFilter->SetSigmaGradient(1.3);
    //houghFilter->SetVariance(2);
    //
    //houghFilter->SetDiscRadiusRatio(0);


    houghFilter->Update();
    AccumulatorImageType::Pointer localAccumulator = houghFilter->GetOutput();
    // Software Guide : EndCodeSnippet
    //  Software Guide : BeginLatex
    //
    //  We can also get the circles as \doxygen{EllipseSpatialObject}. The
    //  \code{GetCircles()} function return a list of those.
    //
    //  Software Guide : EndLatex
    // Software Guide : BeginCodeSnippet
    HoughTransformFilterType::CirclesListType circles;
    circles = houghFilter->GetCircles();
    std::cout << "Found " << circles.size() << " circle(s)." << std::endl;
    // Software Guide : EndCodeSnippet
    //  Software Guide : BeginLatex
    //
    //  We can then allocate an image to draw the resulting circles as binary
    //  objects.
    //
    //  Software Guide : EndLatex
    // Software Guide : BeginCodeSnippet
    typedef  unsigned char                            OutputPixelType;
    typedef  itk::Image< PixelType, 2 > OutputImageType;
    ImageType2D::Pointer  localOutputImage = ImageType2D::New();
    ImageType2D::RegionType region;
    region.SetSize(input->GetOutput()->GetLargestPossibleRegion().GetSize());
    region.SetIndex(input->GetOutput()->GetLargestPossibleRegion().GetIndex());
    localOutputImage->SetRegions(region);
    localOutputImage->SetOrigin(input->GetOutput()->GetOrigin());
    localOutputImage->SetSpacing(input->GetOutput()->GetSpacing());
    localOutputImage->Allocate(true); // initializes buffer to zero
                                      // Software Guide : EndCodeSnippet
                                      //  Software Guide : BeginLatex
                                      //
                                      //  We iterate through the list of circles and we draw them.
                                      //
                                      //  Software Guide : EndLatex
                                      // Software Guide : BeginCodeSnippet
    
    CirclesListType::const_iterator itCircles = circles.begin();
    while (itCircles != circles.end())
    {
        std::cout << "Center: ";
        std::cout << (*itCircles)->GetObjectToParentTransform()->GetOffset()
            << std::endl;
        std::cout << "Radius: " << (*itCircles)->GetRadius()[0] << std::endl;

        if (!keep_only_inside)
        {
            // Software Guide : EndCodeSnippet
            //  Software Guide : BeginLatex
            //
            //  We draw white pixels in the output image to represent each circle.
            //
            //  Software Guide : EndLatex
            // Software Guide : BeginCodeSnippet
            for (double angle = 0;
                angle <= itk::Math::twopi;
                angle += itk::Math::pi / 60.0)
            {
                typedef HoughTransformFilterType::CircleType::TransformType
                    TransformType;
                typedef TransformType::OutputVectorType
                    OffsetType;
                const OffsetType offset =
                    (*itCircles)->GetObjectToParentTransform()->GetOffset();
                localIndex[0] =
                    itk::Math::Round<long int>(offset[0]
                        + (*itCircles)->GetRadius()[0] * std::cos(angle));
                localIndex[1] =
                    itk::Math::Round<long int>(offset[1]
                        + (*itCircles)->GetRadius()[0] * std::sin(angle));
                OutputImageType::RegionType outputRegion =
                    localOutputImage->GetLargestPossibleRegion();
                if (outputRegion.IsInside(localIndex))
                {
                    input->GetOutput()->SetPixel(localIndex, -2);
                }
            }
        }
        else {
            //const ImageType2D * inputImage = input->GetOutput();
            //ImageType2D::RegionType inputRegion = inputImage->GetBufferedRegion();
            //ImageType2D::SizeType size = inputRegion.GetSize();
            //int tolerance = 100;
            //
            //for (unsigned short x = 0; x < size[0]; x++) {
            //    for (unsigned short y = 0; y < size[1]; y++) {
            //        int x_cube = (x - (*itCircles)->GetObjectToParentTransform()->GetOffset()[0])*(x - (*itCircles)->GetObjectToParentTransform()->GetOffset()[0]);
            //        int y_cube = (y - (*itCircles)->GetObjectToParentTransform()->GetOffset()[1])*(y - (*itCircles)->GetObjectToParentTransform()->GetOffset()[1]);
            //        if ((x_cube + y_cube) > ((*itCircles)->GetRadius()[0] + tolerance ) * ((*itCircles)->GetRadius()[0]) + tolerance) {
            //            localIndex[0] = x;
            //            localIndex[1] = y;
            //            localIndex3D[0] = x;
            //            localIndex3D[1] = y;
            //            localIndex3D[2] = actual_slice_number;
            //            input->GetOutput()->SetPixel(localIndex, 0);
            //            original->GetOutput()->SetPixel(localIndex3D, 0);
            //        }
            //    }
            //}

            

        }
        itCircles++;
    }
    
    //#########################################

    /*
    using FilterType2 = itk::FlipImageFilter< ImageType2D >;

    FilterType2::Pointer filter = FilterType2::New();
    // Software Guide : EndCodeSnippet


    //  Software Guide : BeginLatex
    //
    //  The axes to flip are specified in the form of an Array. In this case we
    //  take them from the command line arguments.
    //
    //  \index{itk::FlipImageFilter!Radius}
    //  \index{itk::FlipImageFilter!Neighborhood}
    //
    //  Software Guide : EndLatex

    // Software Guide : BeginCodeSnippet
    using FlipAxesArrayType = FilterType2::FlipAxesArrayType;

    FlipAxesArrayType flipArray;

    flipArray[0] = 1;
    flipArray[1] = 1;

    filter->SetFlipAxes(flipArray);
    // Software Guide : EndCodeSnippet

    //  Software Guide : BeginLatex
    //
    //  The input to the filter can be taken from any other filter, for example
    //  a reader. The output can be passed down the pipeline to other filters,
    //  for example, a writer. Invoking \code{Update()} on any downstream filter
    //  will trigger the execution of the FlipImage filter.
    //
    //  \index{itk::FlipImageFilter!SetInput()}
    //  \index{itk::FlipImageFilter!GetOutput()}
    //
    //  Software Guide : EndLatex


    // Software Guide : BeginCodeSnippet
    //filter->SetInput(localOutputImage);
    filter->SetInput(threshFilter->GetOutput());
    */

    //#########################################
    circles_list.push_back(circles);
    return input;

}


void corrigate_circles_local(std::vector<CirclesListType> &circles_list1) {

    std::vector<CirclesListType> circles_list;

    for (int i = 0; i < circles_list1.size(); i++) {
    
        TransformType::OffsetType Object2ToObject1Offset;
        CirclesListType push_it;
        ElipseType::Pointer elipse = ElipseType::New();
    
        CirclesListType::const_iterator itCircles = circles_list1[i].begin();
        Object2ToObject1Offset[0] = (*itCircles)->GetObjectToParentTransform()->GetOffset()[0];
        Object2ToObject1Offset[1] = (*itCircles)->GetObjectToParentTransform()->GetOffset()[1];
        elipse->GetObjectToParentTransform()->SetOffset(Object2ToObject1Offset);
        elipse->SetRadius((*itCircles)->GetRadius());
        push_it.push_back(elipse);;
        circles_list.push_back(push_it);
    }


    TransformType::OffsetType Object2ToObject1Offset;
    Object2ToObject1Offset[0] = 100;
    Object2ToObject1Offset[1] = 100;


    int filter_size = 15;

    for (int i = 0; i < circles_list.size(); i++) {


        double av_x = 0;
        double av_y = 0;
        double av_rad = 0;

        int sum = 0;
        int j = 0;
        //if (i < std::ceil(filter_size / 2)) {
        //    j = std::ceil(filter_size / 2);
        //}
        //else if (i >= circles_list.size() - std::ceil(filter_size / 2) ) {
        //    j = circles_list.size() - std::ceil(filter_size / 2) - 1;
        //}
        //else {
        //    j = i;
        //}
        
        
        j = i - std::ceil(filter_size / 2);
        int end_of_filter_in_lopp = j + filter_size;
        for (; j < end_of_filter_in_lopp; j++) {
            if (j < 0 || j >= circles_list.size()) {
                continue;
            }
            CirclesListType::const_iterator itCircles = circles_list[j].begin();

            av_x += (*itCircles)->GetObjectToParentTransform()->GetOffset()[0];
            av_y += (*itCircles)->GetObjectToParentTransform()->GetOffset()[1];
            av_rad += (*itCircles)->GetRadius()[0];

            sum++;
        }
        av_x /= sum;
        av_y /= sum;
        av_rad /= sum;

        // AVG
        CirclesListType::const_iterator itCircles = circles_list[i].begin();
        double d = (((*itCircles)->GetObjectToParentTransform()->GetOffset()[0]) - (av_x))*(((*itCircles)->GetObjectToParentTransform()->GetOffset()[0]) - (av_x)) + (((*itCircles)->GetObjectToParentTransform()->GetOffset()[1]) - (av_y))*(((*itCircles)->GetObjectToParentTransform()->GetOffset()[1]) - (av_y));
        if (d < 100000) {
        
        
            CirclesListType::const_iterator first = circles_list1[i].begin();
        
            Object2ToObject1Offset[0] = av_x;
            Object2ToObject1Offset[1] = av_y;
            (*first)->GetObjectToParentTransform()->SetOffset(Object2ToObject1Offset);
            (*first)->SetRadius(av_rad);
        
        
        }

        // MED
        // j = 0;
        // //if (i < std::ceil(filter_size / 2)) {
        // //    j = std::ceil(filter_size / 2);
        // //}
        // //else if (i >= circles_list.size() - std::ceil(filter_size / 2)) {
        // //    j = circles_list.size() - std::ceil(filter_size / 2) - 1;
        // //}
        // //else {
        // //    j = i;
        // //}
        // std::map<int, int> dist_map;
        // j = i - std::ceil(filter_size / 2);
        // end_of_filter_in_lopp = j + filter_size;
        // for (; j < end_of_filter_in_lopp; j++) {
        //     if (j < 0 || j >= circles_list.size()) {
        //         continue;
        //     }
        //     CirclesListType::const_iterator itCircles = circles_list[j].begin();
        // 
        //     int x = (*itCircles)->GetObjectToParentTransform()->GetOffset()[0];
        //     int y = (*itCircles)->GetObjectToParentTransform()->GetOffset()[1];
        //     //double dist = (x - av_x) * (x - av_x) + (y - av_y) * (y - av_y);
        //     double dist = y;
        //     
        //     dist_map[j] = dist;
        // 
        // }
        // 
        // // Declaring the type of Predicate that accepts 2 pairs and return a bool
        // typedef std::function<bool(std::pair<int, int>, std::pair<int, int>)> Comparator;
        // 
        // // Defining a lambda function to compare two pairs. It will compare two pairs using second field
        // Comparator compFunctor =
        //     [](std::pair<int, int> elem1, std::pair<int, int> elem2)
        // {
        //     if (elem1.second == elem2.second) {
        //         return elem1.first < elem2.first;
        //     }
        //     else {
        //         return elem1.second < elem2.second;
        //     }
        // };
        // 
        // // Declaring a set that will store the pairs using above comparision logic
        // std::set<std::pair<int, int>, Comparator> ordered_set(
        //     dist_map.begin(), dist_map.end(), compFunctor);
        // 
        // 
        // 
        // // Iterate over a set using range base for loop
        // // It will display the items in sorted order of values
        // int med_x;
        // int med_y;
        // int med_rad;
        // int k = 0;
        // 
        // for (std::pair<int, int> temp : ordered_set) {
        //    if (k == std::ceil(ordered_set.size() / 2)) {
        //        CirclesListType::const_iterator itCircles = circles_list[temp.first].begin();
        //        med_x = (*itCircles)->GetObjectToParentTransform()->GetOffset()[0];
        //        med_y = (*itCircles)->GetObjectToParentTransform()->GetOffset()[1];
        //        med_rad = (*itCircles)->GetRadius()[0];
        //        dist_map.clear();
        //        break;
        //    }
        //    k++;
        // }
        // 
        // CirclesListType::const_iterator first = circles_list1[i].begin();
        // Object2ToObject1Offset[0] = med_x;
        // Object2ToObject1Offset[1] = med_y;
        // (*first)->GetObjectToParentTransform()->SetOffset(Object2ToObject1Offset);
        // (*first)->SetRadius(med_rad);

    }
    


}

ImageType::Pointer close_on_2d_slices(ImageType::Pointer input) {
    std::cout << "Closing" << std::endl;
    using StructuringElementType = itk::BinaryBallStructuringElement<ImageType::PixelType, ImageType::ImageDimension>;
    StructuringElementType structuringElement;
    structuringElement.SetRadius(1);
    structuringElement.CreateStructuringElement();

    using BinaryMorphologicalClosingImageFilterType =
        itk::BinaryMorphologicalClosingImageFilter<ImageType, ImageType, StructuringElementType>;
    BinaryMorphologicalClosingImageFilterType::Pointer closingFilter = BinaryMorphologicalClosingImageFilterType::New();
    closingFilter->SetForegroundValue(1);
    closingFilter->SetKernel(structuringElement);

    using PasteFilterType = itk::PasteImageFilter<ImageType, ImageType>;
    PasteFilterType::Pointer pasteFilter = PasteFilterType::New();

    const ImageType * inputImage = input;

    ImageType::RegionType inputRegion = inputImage->GetBufferedRegion();
    ImageType::SizeType size = inputRegion.GetSize();
    int height_of_the_image;
    height_of_the_image = size[2];
    using ExtractFilterType = itk::ExtractImageFilter< ImageType, ImageType >;

    pasteFilter->SetDestinationImage(inputImage);
    for (int i = 1; i < height_of_the_image; i++) {

        ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
        extractFilter->SetDirectionCollapseToSubmatrix();

        size[2] = 1; // we extract along z direction

        ImageType::IndexType start = inputRegion.GetIndex();

        start[2] = i;
        ImageType::RegionType desiredRegion;
        desiredRegion.SetSize(size);
        desiredRegion.SetIndex(start);
        extractFilter->SetExtractionRegion(desiredRegion);;
        extractFilter->SetInput(inputImage);

        closingFilter->SetInput(extractFilter->GetOutput());
        pasteFilter->SetSourceImage(closingFilter->GetOutput());

        ImageType::SizeType indexRadius;
        indexRadius[0] = 1; // radius along x
        indexRadius[1] = 1; // radius along y
        indexRadius[2] = 0; // radius along z
        closingFilter->SetRadius(indexRadius);
        closingFilter->UpdateLargestPossibleRegion();
        closingFilter->Update();
        const ImageType * closingImage = closingFilter->GetOutput();
        pasteFilter->SetSourceRegion(closingImage->GetBufferedRegion());
        pasteFilter->SetDestinationIndex(start);

        pasteFilter->Update();
        pasteFilter->SetDestinationImage(pasteFilter->GetOutput());
    }
    std::cout << "Closing end" << std::endl;
    return pasteFilter->GetOutput();
}

template<class MyType>
void normalize_old(MyType input) {
    ImageType::IndexType index3d;
    for (int x = 0; x < input->GetOutput()->GetBufferedRegion().GetSize()[0]; x++) {
        for (int y = 0; y < input->GetOutput()->GetBufferedRegion().GetSize()[1]; y++) {
            for (int z = 0; z < input->GetOutput()->GetBufferedRegion().GetSize()[2]; z++) {
                index3d[0] = x;
                index3d[1] = y;
                index3d[2] = z;

                PixelType pix = input->GetOutput()->GetPixel(index3d);

                if (pix < -1024) {
                    input->GetOutput()->SetPixel(index3d, -1024);
                }
                if (pix > 3071) {
                    input->GetOutput()->SetPixel(index3d, 3071);
                }
            }
        }
    }
}

ScalarImageToHistogramGeneratorType::Pointer normalize(ImageType::Pointer input) {

    ScalarImageToHistogramGeneratorType::Pointer scalarImageToHistogramGenerator = ScalarImageToHistogramGeneratorType::New();

    scalarImageToHistogramGenerator->SetHistogramMin(0);
    scalarImageToHistogramGenerator->SetHistogramMax(1000);
    scalarImageToHistogramGenerator->SetNumberOfBins(50);
    scalarImageToHistogramGenerator->SetAutoHistogramMinimumMaximum(false);
    scalarImageToHistogramGenerator->SetInput(input);
    scalarImageToHistogramGenerator->Compute();

    return scalarImageToHistogramGenerator;
  
}


CalculatorType::OutputType otusHist(ScalarImageToHistogramGeneratorType::Pointer input) {

    //OtsuFilterType::Pointer otsuFilter = OtsuFilterType::New();
    CalculatorType::Pointer otsuFilter = CalculatorType::New();
    //otsuFilter->SetInput(reader_original->GetOutput());
    otsuFilter->SetInputHistogram(input->GetOutput());
    otsuFilter->SetNumberOfThresholds(1);
    otsuFilter->Update(); // To compute threshold

                          //OtsuFilterType::ThresholdVectorType thresholds = otsuFilter->GetThresholds();
    const CalculatorType::OutputType & thresholds = otsuFilter->GetOutput();
    return thresholds;
}

OtsuFilterType::ThresholdVectorType otusImg(ReaderType::Pointer input) {

    OtsuFilterType::Pointer otsuFilter = OtsuFilterType::New();
    otsuFilter->SetInput(input->GetOutput());
    otsuFilter->SetNumberOfThresholds(1);
    otsuFilter->Update(); // To compute threshold

    OtsuFilterType::ThresholdVectorType thresholds = otsuFilter->GetThresholds();
    return thresholds;
}

//std::string files[] = { "test\\0verse112.nii.gz", "test\\1AnyScan.nii.gz", "test\\2LA_LD_02.nii.gz", "test\\3LA_1.nii.gz", "test\\4LUNG_PET.nii.gz", "test\\5LungRT.nii.gz", "test\\6LA_Diagn_1.nii.gz", "test\\7LA_Diagn_3.nii.gz", "test\\8LA_LD_01.nii.gz", "test\\9LungPhilips.nii.gz" };
//std::string files[] = { "test\\2LA_LD_02.nii.gz" };
//std::string files[] = { "test\\0verse112.nii.gz" };


//std::string files[] = { "bigtest\\lung_08.nii.gz" };

int main(int argc, char ** argv)
{

    std::vector<std::string> files;
    for (int i = 1; i <= 60; i++) {
        if (i > 9) {
            files.push_back("bigtest\\lung_" + std::to_string(i) + ".nii.gz");
    
        }
        else {
            files.push_back("bigtest\\lung_0" + std::to_string(i) + ".nii.gz");
        }
    }

   /* if (argv[1] != nullptr) {
        files = { argv[1] };
    }*/

    int number_of_test_file = 1;
    for (std::string file : files) {
        if ("test\\8LA_LD_01.nii.gz" == file) {
            number_of_test_file++;
            continue;
        }
        using WriterTypeT = itk::ImageFileWriter< ImageType >;
        WriterTypeT::Pointer writerT = WriterTypeT::New();
        //// Verify the number of parameters in the command line
        //if (argc <= 3)
        //{
        //    std::cerr << "Usage: " << std::endl;
        //    std::cerr << argv[0] << " input3DImageFile  output3DImageFile " << std::endl;
        //    std::cerr << " sliceNumber " << std::endl;
        //    return EXIT_FAILURE;
        //}

        //const char * original_image_file_name = argv[1];
        std::string original_image_file_name = file;

        // Here we recover the file names from the command line arguments
        /*const char * inputFilename = argv[1];
    

        ReaderType::Pointer reader = ReaderType::New();
        reader->SetFileName(inputFilename);
        try
        {
            reader->Update();
        }
        catch (itk::ExceptionObject & err)
        {
            std::cerr << "ExceptionObject caught !" << std::endl;
            std::cerr << err << std::endl;
            return EXIT_FAILURE;
        }*/

        ReaderType::Pointer reader_original = ReaderType::New();
        reader_original->SetFileName(original_image_file_name);
        try
        {
            reader_original->Update();
        }
        catch (itk::ExceptionObject & err)
        {
            std::cerr << "ExceptionObject caught !" << std::endl;
            std::cerr << err << std::endl;
            return EXIT_FAILURE;
        }

        ReaderType::Pointer reader_original_1 = ReaderType::New();
        reader_original_1->SetFileName(original_image_file_name);
        try
        {
            reader_original_1->Update();
        }
        catch (itk::ExceptionObject & err)
        {
            std::cerr << "ExceptionObject caught !" << std::endl;
            std::cerr << err << std::endl;
            return EXIT_FAILURE;
        }

        ReaderType::Pointer reader_original_2 = ReaderType::New();
        reader_original_2->SetFileName(original_image_file_name);
        try
        {
            reader_original_2->Update();
        }
        catch (itk::ExceptionObject & err)
        {
            std::cerr << "ExceptionObject caught !" << std::endl;
            std::cerr << err << std::endl;
            return EXIT_FAILURE;
        }

        ScalarImageToHistogramGeneratorType::Pointer a = normalize(reader_original->GetOutput());
        ScalarImageToHistogramGeneratorType::Pointer b = normalize(reader_original_1->GetOutput());
        ScalarImageToHistogramGeneratorType::Pointer c = normalize(reader_original_2->GetOutput());

        // set up the extraction region [one slice]
    
        //const ImageType * inputImage = reader->GetOutput();

        std::cout << "Thresholding" << std::endl;

        //OtsuFilterType::Pointer otsuFilter = OtsuFilterType::New();
        //otsuFilter->SetInput(reader_original->GetOutput());
        //otsuFilter->SetNumberOfThresholds(4);
        //otsuFilter->Update(); // To compute threshold

        //OtsuFilterType::ThresholdVectorType thresholds = otsuFilter->GetThresholds();
        const CalculatorType::OutputType & thresholds = otusHist(a);
        for (unsigned int i = 0; i < thresholds.size(); i++)
        {
            std::cout << thresholds[i] << std::endl;
        }

        using FilterType = itk::BinaryThresholdImageFilter< ImageType, ImageType >;
        FilterType::Pointer threshFilter = FilterType::New();
        threshFilter->SetInput(reader_original->GetOutput());
        //if (thresholds[3] > 1000) {
        //    threshFilter->SetLowerThreshold(thresholds[2]);
        //}
        //else {
        //    threshFilter->SetLowerThreshold(thresholds[3]);
        //}
        threshFilter->SetLowerThreshold(thresholds[0]);
        threshFilter->SetUpperThreshold(10000);
        threshFilter->SetOutsideValue(0);
        threshFilter->SetInsideValue(1);
        threshFilter->Update();

        // CLOSING
        const ImageType::Pointer inputImage = close_on_2d_slices(threshFilter->GetOutput());
        //const ImageType * inputImage = threshFilter->GetOutput();

        ImageType::RegionType inputRegion = inputImage->GetBufferedRegion();
        ImageType::SizeType size = inputRegion.GetSize();
        int height_of_the_image;
        height_of_the_image = size[2];
        using ExtractFilterType = itk::ExtractImageFilter< ImageType, ImageType2D >;
        std::list<ExtractFilterType::Pointer> slice_list;
        std::vector<CirclesListType> circles_list;
        for (int i = 0; i < height_of_the_image; i++) {

            ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
            extractFilter->SetDirectionCollapseToSubmatrix();


            size[2] = 0; // we extract along z direction

            ImageType::IndexType start = inputRegion.GetIndex();

            //const unsigned int sliceNumber = std::stoi( argv[3] );

            start[2] = i;
            ImageType::RegionType desiredRegion;
            desiredRegion.SetSize(size);
            desiredRegion.SetIndex(start);

            extractFilter->SetExtractionRegion(desiredRegion);

            extractFilter->SetInput(inputImage); 

            slice_list.push_back(do_hough_on_image<ExtractFilterType::Pointer, ReaderType::Pointer>(extractFilter, reader_original, true, i, circles_list));

        }

        bool zero_c_flag = false;
        auto it = circles_list.begin();
        while (it != circles_list.end())
        {
            // CirclesListType::const_iterator itCircles = (*it).begin();
            int s = (*it).size();
            // std::cout << s << std::endl;
            if (s < 1) {
                zero_c_flag = true;
                break;
            }
            it++;
        }
        if (zero_c_flag) {
            number_of_test_file++;
            continue;
        }

        draw_circle<ReaderType::Pointer>(reader_original, circles_list);
        std::string name_of_the_file = std::to_string(number_of_test_file) + "_original_circles_0_corrig.nii.gz";

        writerT->SetFileName(name_of_the_file);
        writerT->SetInput(reader_original->GetOutput());
        try
        {
            writerT->Update();
        }
        catch (itk::ExceptionObject & err)
        {
            std::cerr << "ExceptionObject caught !" << std::endl;
            std::cerr << err << std::endl;
            return EXIT_FAILURE;
        }

        corrigate_circles(circles_list);
        draw_circle<ReaderType::Pointer>(reader_original_1, circles_list);

        name_of_the_file = std::to_string(number_of_test_file) + "_original_circles_1_corrig.nii.gz";

        writerT->SetFileName(name_of_the_file);
        writerT->SetInput(reader_original_1->GetOutput());
        try
        {
            writerT->Update();
        }
        catch (itk::ExceptionObject & err)
        {
            std::cerr << "ExceptionObject caught !" << std::endl;
            std::cerr << err << std::endl;
            return EXIT_FAILURE;
        }


        corrigate_circles_local(circles_list);
        draw_circle<ReaderType::Pointer>(reader_original_2, circles_list);
        //draw_circle<FilterType::Pointer>(threshFilter, circles_list);

        name_of_the_file = std::to_string(number_of_test_file) + "_original_circles_2_corrig.nii.gz";
        
        writerT->SetFileName(name_of_the_file);
        writerT->SetInput(reader_original_2->GetOutput());
        try
        {
            writerT->Update();
        }
        catch (itk::ExceptionObject & err)
        {
            std::cerr << "ExceptionObject caught !" << std::endl;
            std::cerr << err << std::endl;
            return EXIT_FAILURE;
        }


        //name_of_the_file = std::to_string(number_of_test_file) + "_otsu_circles.nii.gz";
        //writerT->SetFileName(name_of_the_file);
        //writerT->SetInput(threshFilter->GetOutput());
        //try
        //{
        //    writerT->Update();
        //}
        //catch (itk::ExceptionObject & err)
        //{
        //    std::cerr << "ExceptionObject caught !" << std::endl;
        //    std::cerr << err << std::endl;
        //    return EXIT_FAILURE;
        //}
        //
        //// OTSU ON CIRCLES
        //otsuFilter->SetInput(reader_original->GetOutput());
        //otsuFilter->SetNumberOfThresholds(2);
        //otsuFilter->Update(); // To compute threshold
        //
        //thresholds = otsuFilter->GetThresholds();
        //for (unsigned int i = 0; i < thresholds.size(); i++)
        //{
        //    std::cout << thresholds[i] << std::endl;
        //}
        //
        //
        //threshFilter->SetInput(reader_original->GetOutput());
        //threshFilter->SetLowerThreshold(thresholds[1]);
        //threshFilter->SetUpperThreshold(10000);
        //threshFilter->SetOutsideValue(0);
        //threshFilter->SetInsideValue(1);
        //threshFilter->Update();
        //
        //name_of_the_file = std::to_string(number_of_test_file) + "_thresh_on_circles.nii.gz";
        //writerT->SetFileName(name_of_the_file);
        //writerT->SetInput(threshFilter->GetOutput());
        //try
        //{
        //    writerT->Update();
        //}
        //catch (itk::ExceptionObject & err)
        //{
        //    std::cerr << "ExceptionObject caught !" << std::endl;
        //    std::cerr << err << std::endl;
        //    return EXIT_FAILURE;
        //}


        number_of_test_file++;
    }

    return EXIT_SUCCESS;
}