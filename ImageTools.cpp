// ImageTools contain some small tools for medical image processing
// Write by Noah Jing Li based on Boost library imitated FISH library.

#include "stdafx.h"
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

//using namespace UNIBE::IPMI;

#define TYPE short //Noah comment: probility map use "double" while anatomical/label image use "short"

void print_help()
{
	//std::cerr << "[ dataTypeD | int, uint, short, ushort, float, double, char, uchar ]\n" << std::endl;

	std::cerr << "[ nearest upsample image ]" << std::endl;
	std::cerr << "example: ImageTools --nearest_upsample --ref referenceImage --in inputImage --out outputImage --spacing value\n" << std::endl;	

	std::cerr << "[ linear upsample image ]" << std::endl;
	std::cerr << "example: ImageTools --linear_upsample --ref referenceImage --in inputImage --out outputImage --spacing value\n" << std::endl;

	std::cerr << "[ upsample image by itself ]" << std::endl;
	std::cerr << "example: ImageTools --upsample_self --in inputImage --out outputImage --spacing_file spacing_file\n" << std::endl;

	std::cerr << "[ label upsample image by itself ]" << std::endl;
	std::cerr << "example: ImageTools --label_upsample_self --in inputImage --out outputImage --spacing_file spacing_file\n" << std::endl;

	std::cerr << "[ resample image ]" << std::endl;
	std::cerr << "example: ImageTools --resample --in inputImage --out outputImage --spacing_file spacing_file\n" << std::endl;

	std::cerr << "[ downsample image ]" << std::endl;
	std::cerr << "example: ImageTools --downsample --in inputImage --spacing value --out outputImage\n" << std::endl;

	std::cerr << "[ generate mask image ]" << std::endl;
	std::cerr << "example: ImageTools --mask --in inputImage --threshold value --out outputMaskImage\n" << std::endl;

	std::cerr << "[ convert point from world coord to voxel coord ]" << std::endl;
	std::cerr << "In the pointfile, point should be split by Tab." << std::endl;
	std::cerr << "example: ImageTools --world2voxel --in inputImage --landmark pointfile\n" << std::endl;
	
	std::cerr << "[ convert point from voxel coord to world coord ]" << std::endl;
	std::cerr << "In the pointfile, point should be split by Tab." << std::endl;
	std::cerr << "example: ImageTools --voxel2world --in inputImage --landmark pointfile\n" << std::endl;
	
	std::cerr << "[ histogram matching ]" << std::endl;
	std::cerr << "example: ImageTools --histogram_match --ref referenceImage --in inputImage --out outputImage\n" << std::endl;

	std::cerr << "[ generate the segmentation ]" << std::endl;
	std::cerr << "example: ImageTools --generate_segment --in probabilityMap --threshold value --out segmentationImage\n" << std::endl;

	std::cerr << "[ extract largest connected component from label image ]" << std::endl;
	std::cerr << "example: ImageTools --largest_component --in inputImage --out outputImage\n" << std::endl;

	std::cerr << "[ fill the holes in the label image ]" << std::endl;
	std::cerr << "example: ImageTools --fill_hole --in inputImage --out outputImage --radius value\n" << std::endl;
	
	std::cerr << "[ (under development) correct the image to origin(0,0,0), direction(1,1,1), isotropy spacing ]" << std::endl;
	std::cerr << "example: ImageTools --correct --ref referenceImage --in inputImage --out outputImage --spacing value\n" << std::endl;

	std::cerr << "[ dilate binary image ]" << std::endl;
	std::cerr << "example: ImageTools --dilate --in inputImage --radius value --out outputImage\n" << std::endl;

	std::cerr << "[ erode binary image ]" << std::endl;
	std::cerr << "example: ImageTools --erode --in inputImage --radius value --out outputImage\n" << std::endl;

	std::cerr << "[ convert data type to 'nii.gz' ]" << std::endl;
	std::cerr << "example: ImageTools --convert --in inputImage --out outputImage\n" << std::endl;

	std::cerr << "[ count voxel number ]" << std::endl;
	std::cerr << "example: ImageTools --count --in inputImage\n" << std::endl;

	std::cerr << "[ crop volumn ]" << std::endl;
	std::cerr << "example: ImageTools --crop --in inputImage --start_pt filename --end_pt filename --pad value mm --out outputImage\n" << std::endl;

	std::cerr << "[ crop label volume for femur based on lesser trochanter landmark ]" << std::endl;
	std::cerr << "example: ImageTools --crop_femur --in inputImage --landmark pointfile --out outputImage\n" << std::endl;

	std::cerr << "[ get bounding box ]" << std::endl;
	std::cerr << "example: ImageTools --boundingbox --in inputImage --start_pt filename --end_pt filename\n" << std::endl;

	std::cerr << "[ image minus specified value ]" << std::endl;
	std::cerr << "example: ImageTools --minus --in inputImage --minus_value value --out outputImage\n" << std::endl;

	std::cerr << "[ put cropped image back to original image ]" << std::endl;
	std::cerr << "example: ImageTools --putback --ref referenceImage --in inputImage --start_pt filename --end_pt filename --out outputImage\n" << std::endl;

	std::cerr << "[ image label transform ]" << std::endl;
	std::cerr << "example: ImageTools --label_transform --in inputImage --out outputImage\n" << std::endl;

	std::cerr << "[ adaptive histogram equalization ]" << std::endl;
	std::cerr << "example: ImageTools --adaptive_histogram_equalization --in inputImage --radius value --alpha value --beta value --out outputImage\n" << std::endl;

	std::cerr << "[ rescale intensity ]" << std::endl;
	std::cerr << "example: ImageTools --rescale_intensity --in inputImage --min value --max value --out outputImage\n" << std::endl;

	std::cerr << "[ (under development) orient image to specified direction ]" << std::endl;
	std::cerr << "example: ImageTools --orient --in inputImage --out outputImage\n" << std::endl;

	std::cerr << "[ judge two image is same or not ]" << std::endl;
	std::cerr << "example: ImageTools --judge --ref referenceImage --in inputImage \n" << std::endl;

	std::cerr << "[ cut off the histogram ]" << std::endl;
	std::cerr << "example: ImageTools --cutoff --in inputImage --out outputImage --ratio value \n" << std::endl;

	std::cerr << "[ extract left-top and right-bottom landmarks for femur ROI based on bounding box ]" << std::endl;
	std::cerr << "example: ImageTools --extract_landmark_femur --in inputLabelImage --start_pt filename --end_pt filename  --out outputMaskImage \n" << std::endl;
}

//void CreateImage(ImageType::Pointer image)
//{
//	// Create an image with 2 connected components
//	ImageType::RegionType region;
//	ImageType::IndexType start;
//	start[0] = 0;
//	start[1] = 0;
//
//	ImageType::SizeType size;
//	unsigned int NumRows = 200;
//	unsigned int NumCols = 300;
//	size[0] = NumRows;
//	size[1] = NumCols;
//
//	region.SetSize(size);
//	region.SetIndex(start);
//
//	image->SetRegions(region);
//	image->Allocate();
//	image->FillBuffer(itk::NumericTraits<ImageType::PixelType>::Zero);
//
//	// Make a rectangle, centered at (100,150) with sides 160 & 240
//	// This provides a 20 x 30 border around the square for the crop filter to remove
//	for (unsigned int r = 20; r < 180; r++)
//	{
//		for (unsigned int c = 30; c < 270; c++)
//		{
//			ImageType::IndexType pixelIndex;
//			pixelIndex[0] = r;
//			pixelIndex[1] = c;
//
//			image->SetPixel(pixelIndex, 200);
//		}
//	}
//}

enum LoadPixelType {
	UNKNOWN, _UCHAR, _CHAR, _USHORT, _SHORT, _UINT, _INT, _ULONG, _LONG, _FLOAT, _DOUBLE
};

/** @brief read pixel type */
LoadPixelType ReadPixelType(const char* path) {

	typedef itk::ImageIOBase::IOComponentType ScalarPixelType;
	itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(path, itk::ImageIOFactory::ReadMode);
	if (imageIO.IsNull()) {
		std::cerr << "image file does not exist: " << path << std::endl; exit(-1);
	}

	imageIO->SetFileName(path);
	imageIO->ReadImageInformation();

	if (imageIO->GetPixelType() != itk::ImageIOBase::SCALAR) {
		std::cerr << "the pixel type is not a scalar" << std::endl; exit(-1);
	}

	const ScalarPixelType pixelType = imageIO->GetComponentType();
	switch (pixelType)
	{
	case itk::ImageIOBase::UCHAR: return _UCHAR; break;
	case itk::ImageIOBase::CHAR: return _CHAR; break;
	case itk::ImageIOBase::USHORT: return _USHORT; break;
	case itk::ImageIOBase::SHORT: return _SHORT; break;
	case itk::ImageIOBase::UINT: return _UINT; break;
	case itk::ImageIOBase::INT: return _INT; break;
	case itk::ImageIOBase::ULONG: return _ULONG; break;
	case itk::ImageIOBase::LONG: return _LONG; break;
	case itk::ImageIOBase::FLOAT: return _FLOAT; break;
	case itk::ImageIOBase::DOUBLE: return _DOUBLE; break;
	default:
		return UNKNOWN;
		break;
	}
}


template <typename T>
T Get(boost::program_options::variables_map& vm, const char* name)
{
	T val;
	if (vm.count(name)) {
		val = vm[name].as<T>();
	}
	else {
		std::cerr << name << " option missing" << std::endl;
		exit(-1);
	}
	return val;
}

//under development
template <typename InputType>
itk::Image<InputType, 3>* ReadImage(const char* path) {
	typedef itk::Image<InputType, 3>					ImageType;
	typedef itk::ImageFileReader< ImageType >			ReaderType;

	typename ReaderType::Pointer reader = ReaderType::New();

	reader->SetFileName(path);
	try {
		reader->Update();
	}
	catch (const itk::ExceptionObject& exp) {
		std::cerr << exp.GetDescription() << std::endl;
		return nullptr;
	}

	return reader->GetOutput();
}

//under development
template <typename InputType, typename OutputType>
bool WriteImage(typename itk::Image<InputType, 3>::Pointer image, const char* path) {
	
	//typedef typename itk::Image<OutputType, 3>								InternalImageType;
	////typedef typename itk::CastImageFilter<InputType, OutputType>				ImageCasterType;
	////
	////typename ImageCasterType::Pointer imageCaster = ImageCasterType::New();
	////imageCaster->SetInput(image);

	//typedef typename itk::ImageFileWriter<InternalImageType>					FileWriterType;
	//typename FileWriterType::Pointer writer = FileWriterType::New();

	//writer->SetFileName(path);
	//writer->SetInput(image);//error
	//try {
	//	writer->Update();
	//}
	//catch (const itk::ExceptionObject& exp) {
	//	std::cerr << exp.GetDescription() << std::endl;
	//	return false;
	//}

	return true;

}

///////////////////////////////////////////////////////////////////////////////
//Goal: nearest upsample image
//Step: 1. load the reference image
//		2. load the upsample image candidate
//		3. BSpline/Linear upsample 
//		4. save the upsampled image
///////////////////////////////////////////////////////////////////////////////
template<typename T>
bool Nearest_UpSample(const std::string referenceImage, const std::string inputImageFile, const std::string outputImageFile, double spacing_value) {
	
	//Step: 1. load the reference image
	typedef T												PixelType;
	const unsigned int Dimension = 3;

	typedef itk::Image<PixelType, Dimension>				ImageType;
	typedef itk::ImageFileReader<ImageType>					ReaderType;

	//typedef T 												T_WritePixel;
	typedef itk::Image<PixelType, Dimension>				WriterImageType;
	typedef itk::ImageFileWriter<WriterImageType>				WriterType;

	ReaderType::Pointer reader_refer = ReaderType::New();

	const char* imagefile_refer = referenceImage.c_str();
	if (!imagefile_refer) {
		std::cout << "Load input file error!" << std::endl;
		return EXIT_FAILURE;
	}

	reader_refer->SetFileName(imagefile_refer);
	reader_refer->Update();

	std::cout << "Read " << imagefile_refer << " DONE!" << std::endl;


	//2. load the upsample image candidate
	ReaderType::Pointer reader_candidate = ReaderType::New();

	const char* imagefile_candidate = inputImageFile.c_str();
	if (!imagefile_candidate) {
		std::cout << "Load input file error!" << std::endl;
		return EXIT_FAILURE;
	}

	reader_candidate->SetFileName(imagefile_candidate);
	reader_candidate->Update();

	std::cout << "Read " << imagefile_candidate << " DONE!" << std::endl;

	ImageType::SpacingType spacing;
	spacing[0] = spacing_value;
	spacing[1] = spacing_value;
	spacing[2] = spacing_value;

	// 3. BSpline upsample
	typedef itk::IdentityTransform<double, Dimension>								TransformType;

	//typedef itk::BSplineInterpolateImageFunction<ImageType, double, TYPE>			InterpolatorType;
	//typedef itk::LinearInterpolateImageFunction< ImageType, double >				InterpolatorType;
	typedef itk::NearestNeighborInterpolateImageFunction< ImageType, double >		InterpolatorType;
	
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	//interpolator->SetSplineOrder(3);// only for Bspline interpolate

	typedef itk::ResampleImageFilter<ImageType, ImageType>						ResampleFilterType;

	TransformType::Pointer transform = TransformType::New();
	transform->SetIdentity();

	ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();
	resampleFilter->SetTransform(transform);
	resampleFilter->SetInterpolator(interpolator);
	resampleFilter->SetOutputOrigin(reader_refer->GetOutput()->GetOrigin());
	resampleFilter->SetOutputSpacing(spacing);
	resampleFilter->SetSize(reader_refer->GetOutput()->GetLargestPossibleRegion().GetSize());
	resampleFilter->SetOutputDirection(reader_refer->GetOutput()->GetDirection());
	resampleFilter->SetDefaultPixelValue(0);
	resampleFilter->SetInput(reader_candidate->GetOutput());

	// 4. save the upsampled image
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(outputImageFile.c_str());
	writer->SetInput(resampleFilter->GetOutput());
	writer->Update();

	return EXIT_SUCCESS;
}
//bool Nearest_UpSample(const std::string referenceImage, const std::string inputImageFile, const std::string outputImageFile, double spacing_value) {
//	//For experiment, not the standard version.
//	//Step: 1. load the reference image
//	typedef TYPE												PixelType;
//	const unsigned int Dimension = 3;
//
//	typedef itk::Image<PixelType, Dimension>				ImageType;
//	typedef itk::ImageFileReader<ImageType>					ReaderType;
//
//	typedef TYPE T_WritePixel;
//	typedef itk::Image<T_WritePixel, Dimension>				WriterImageType;
//	typedef itk::ImageFileWriter<WriterImageType>				WriterType;
//
//	ReaderType::Pointer reader_refer = ReaderType::New();
//
//	const char* imagefile_refer = referenceImage.c_str();
//	if (!imagefile_refer) {
//		std::cout << "Load input file error!" << std::endl;
//		return EXIT_FAILURE;
//	}
//
//	reader_refer->SetFileName(imagefile_refer);
//	reader_refer->Update();
//
//	std::cout << "Read " << imagefile_refer << " DONE!" << std::endl;
//
//
//	//2. load the upsample image candidate
//	ReaderType::Pointer reader_candidate = ReaderType::New();
//
//	const char* imagefile_candidate = inputImageFile.c_str();
//	if (!imagefile_candidate) {
//		std::cout << "Load input file error!" << std::endl;
//		return EXIT_FAILURE;
//	}
//
//	reader_candidate->SetFileName(imagefile_candidate);
//	reader_candidate->Update();
//
//	std::cout << "Read " << imagefile_candidate << " DONE!" << std::endl;
//
//	ImageType::SpacingType spacing;
//	spacing[0] = spacing_value;
//	spacing[1] = spacing_value;
//	spacing[2] = spacing_value;
//
//	// 3. BSpline upsample
//	typedef itk::IdentityTransform<double, Dimension>								TransformType;
//
//	//typedef itk::BSplineInterpolateImageFunction<ImageType, double, TYPE>			InterpolatorType;
//	//typedef itk::LinearInterpolateImageFunction< ImageType, double >				InterpolatorType;
//	typedef itk::NearestNeighborInterpolateImageFunction< ImageType, TYPE >		InterpolatorType;
//
//	InterpolatorType::Pointer interpolator = InterpolatorType::New();
//	//interpolator->SetSplineOrder(3);// only for Bspline interpolate
//
//	typedef itk::ResampleImageFilter<ImageType, ImageType>						ResampleFilterType;
//
//	TransformType::Pointer transform = TransformType::New();
//	transform->SetIdentity();
//
//	ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();
//	resampleFilter->SetTransform(transform);
//	resampleFilter->SetInterpolator(interpolator);
//	resampleFilter->SetOutputOrigin(reader_refer->GetOutput()->GetOrigin());
//	resampleFilter->SetOutputSpacing(reader_candidate->GetOutput()->GetSpacing());
//	resampleFilter->SetSize(reader_candidate->GetOutput()->GetLargestPossibleRegion().GetSize());
//	resampleFilter->SetOutputDirection(reader_refer->GetOutput()->GetDirection());
//	resampleFilter->SetDefaultPixelValue(0);
//	resampleFilter->SetInput(reader_candidate->GetOutput());
//
//	// 4. save the upsampled image
//	WriterType::Pointer writer = WriterType::New();
//	writer->SetFileName(outputImageFile.c_str());
//	writer->SetInput(resampleFilter->GetOutput());
//	writer->Update();
//
//	return EXIT_SUCCESS;
//}


///////////////////////////////////////////////////////////////////////////////
//Goal: linear upsample image
//Step: 1. load the reference image
//		2. load the upsample image candidate
//		3. BSpline/Linear upsample 
//		4. save the upsampled image
///////////////////////////////////////////////////////////////////////////////
template <typename T>
bool Linear_UpSample(const std::string referenceImage, const std::string inputImageFile, const std::string outputImageFile, double spacing_value) {

	//Step: 1. load the reference image
	typedef T												PixelType;
	const unsigned int Dimension = 3;

	typedef itk::Image<PixelType, Dimension>				ImageType;
	typedef itk::ImageFileReader<ImageType>					ReaderType;

	typedef T 												T_WritePixel;
	typedef itk::Image<T_WritePixel, Dimension>				WriterImageType;
	typedef itk::ImageFileWriter<WriterImageType>				WriterType;

	ReaderType::Pointer reader_refer = ReaderType::New();

	const char* imagefile_refer = referenceImage.c_str();
	if (!imagefile_refer) {
		std::cout << "Load input file error!" << std::endl;
		return EXIT_FAILURE;
	}

	reader_refer->SetFileName(imagefile_refer);
	reader_refer->Update();

	std::cout << "Read " << imagefile_refer << " DONE!" << std::endl;


	//2. load the upsample image candidate
	ReaderType::Pointer reader_candidate = ReaderType::New();

	const char* imagefile_candidate = inputImageFile.c_str();
	if (!imagefile_candidate) {
		std::cout << "Load input file error!" << std::endl;
		return EXIT_FAILURE;
	}

	reader_candidate->SetFileName(imagefile_candidate);
	reader_candidate->Update();

	std::cout << "Read " << imagefile_candidate << " DONE!" << std::endl;

	ImageType::SpacingType spacing;
	spacing[0] = spacing_value;
	spacing[1] = spacing_value;
	spacing[2] = spacing_value;

	// 3. BSpline upsample
	typedef itk::IdentityTransform<double, Dimension>								TransformType;

	//typedef itk::BSplineInterpolateImageFunction<ImageType, double, TYPE>			InterpolatorType;
	typedef itk::LinearInterpolateImageFunction< ImageType, double >				InterpolatorType;
	//typedef itk::NearestNeighborInterpolateImageFunction< ImageType, TYPE >		InterpolatorType;

	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	//interpolator->SetSplineOrder(3);// only for Bspline interpolate

	typedef itk::ResampleImageFilter<ImageType, ImageType>						ResampleFilterType;

	TransformType::Pointer transform = TransformType::New();
	transform->SetIdentity();

	ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();
	resampleFilter->SetTransform(transform);
	resampleFilter->SetInterpolator(interpolator);
	resampleFilter->SetOutputOrigin(reader_refer->GetOutput()->GetOrigin());
	resampleFilter->SetOutputSpacing(spacing);
	resampleFilter->SetSize(reader_refer->GetOutput()->GetLargestPossibleRegion().GetSize());
	resampleFilter->SetOutputDirection(reader_refer->GetOutput()->GetDirection());
	resampleFilter->SetDefaultPixelValue(0);
	resampleFilter->SetInput(reader_candidate->GetOutput());

	// 4. save the upsampled image
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(outputImageFile.c_str());
	writer->SetInput(resampleFilter->GetOutput());
	writer->Update();

	return EXIT_SUCCESS;
}
//bool Linear_UpSample(const std::string referenceImage, const std::string inputImageFile, const std::string outputImageFile, double spacing_value) {
//	// Only for experiment, not the standard version.
//	//Step: 1. load the reference image
//	typedef TYPE												PixelType;
//	const unsigned int Dimension = 3;
//
//	typedef itk::Image<PixelType, Dimension>				ImageType;
//	typedef itk::ImageFileReader<ImageType>					ReaderType;
//
//	typedef TYPE T_WritePixel;
//	typedef itk::Image<T_WritePixel, Dimension>				WriterImageType;
//	typedef itk::ImageFileWriter<WriterImageType>				WriterType;
//
//	ReaderType::Pointer reader_refer = ReaderType::New();
//
//	const char* imagefile_refer = referenceImage.c_str();
//	if (!imagefile_refer) {
//		std::cout << "Load input file error!" << std::endl;
//		return EXIT_FAILURE;
//	}
//
//	reader_refer->SetFileName(imagefile_refer);
//	reader_refer->Update();
//
//	std::cout << "Read " << imagefile_refer << " DONE!" << std::endl;
//
//
//	//2. load the upsample image candidate
//	ReaderType::Pointer reader_candidate = ReaderType::New();
//
//	const char* imagefile_candidate = inputImageFile.c_str();
//	if (!imagefile_candidate) {
//		std::cout << "Load input file error!" << std::endl;
//		return EXIT_FAILURE;
//	}
//
//	reader_candidate->SetFileName(imagefile_candidate);
//	reader_candidate->Update();
//
//	std::cout << "Read " << imagefile_candidate << " DONE!" << std::endl;
//
//	ImageType::SpacingType spacing;
//	spacing[0] = spacing_value;
//	spacing[1] = spacing_value;
//	spacing[2] = spacing_value;
//
//	// 3. BSpline upsample
//	typedef itk::IdentityTransform<double, Dimension>								TransformType;
//
//	//typedef itk::BSplineInterpolateImageFunction<ImageType, double, TYPE>			InterpolatorType;
//	typedef itk::LinearInterpolateImageFunction< ImageType, double >				InterpolatorType;
//	//typedef itk::NearestNeighborInterpolateImageFunction< ImageType, TYPE >		InterpolatorType;
//
//	InterpolatorType::Pointer interpolator = InterpolatorType::New();
//	//interpolator->SetSplineOrder(3);// only for Bspline interpolate
//
//	typedef itk::ResampleImageFilter<ImageType, ImageType>						ResampleFilterType;
//
//	TransformType::Pointer transform = TransformType::New();
//	transform->SetIdentity();
//
//	ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();
//	resampleFilter->SetTransform(transform);
//	resampleFilter->SetInterpolator(interpolator);
//	resampleFilter->SetOutputOrigin(reader_refer->GetOutput()->GetOrigin());
//	resampleFilter->SetOutputSpacing(reader_candidate->GetOutput()->GetSpacing());
//	resampleFilter->SetSize(reader_candidate->GetOutput()->GetLargestPossibleRegion().GetSize());
//	resampleFilter->SetOutputDirection(reader_refer->GetOutput()->GetDirection());
//	resampleFilter->SetDefaultPixelValue(0);
//	resampleFilter->SetInput(reader_candidate->GetOutput());
//
//	// 4. save the upsampled image
//	WriterType::Pointer writer = WriterType::New();
//	writer->SetFileName(outputImageFile.c_str());
//	writer->SetInput(resampleFilter->GetOutput());
//	writer->Update();
//
//	return EXIT_SUCCESS;
//}

///////////////////////////////////////////////////////////////////////////////
//Goal: upsample image by itself
//Step: 1. load the upsample image candidate
//		2. compute the output size
//		3. BSpline/Linear upsample 
//		4. save the upsampled image
///////////////////////////////////////////////////////////////////////////////
template <typename T>
bool UpSample_Self(const std::string input_path, const std::string output_path, const std::string spacing_path) {
	
	//Step: 1. load the upsample image candidate
	typedef T												PixelType;
	const unsigned int Dimension = 3;

	typedef itk::Image<PixelType, Dimension>				ImageType;
	typedef itk::ImageFileReader<ImageType>					ReaderType;

	typedef itk::Image<PixelType, Dimension>				WriterImageType;
	typedef itk::ImageFileWriter<WriterImageType>				WriterType;

	ReaderType::Pointer reader_candidate = ReaderType::New();

	const char* imagefile_candidate = input_path.c_str();
	if (!imagefile_candidate) {
		std::cout << "Load input file error!" << std::endl;
		return EXIT_FAILURE;
	}

	reader_candidate->SetFileName(imagefile_candidate);
	reader_candidate->Update();

	std::cout << "Read " << imagefile_candidate << " DONE!" << std::endl;

	//2. compute the output size
	ImageType::Pointer image = reader_candidate->GetOutput();
	ImageType::RegionType region = image->GetLargestPossibleRegion();
	const ImageType::SizeType size = region.GetSize();
	std::cout << "Input size: " << size << std::endl;

	const ImageType::SpacingType& sp = image->GetSpacing();

	std::cout << "Original Spacing = ";
	std::cout << sp[0] << ", " << sp[1] << ", " << sp[2] << std::endl;

	ImageType::PointType physicalSize;
	for (int i = 0; i < 3; ++i) {
		physicalSize[i] = size[i] * sp[i];
	}

	/*ImageType::SpacingType spacing;
	spacing[0] = spacing_value;
	spacing[1] = spacing_value;
	spacing[2] = spacing_value;*/
	const char* spacingfile = spacing_path.c_str();

	std::ifstream infile;
	infile.open(spacingfile, std::ifstream::in);

	if (!infile) {
		std::cerr << "Load input landmark error!" << std::endl;
		return EXIT_FAILURE;
	}

	std::string stream_spacing;
	getline(infile, stream_spacing);

	infile.close();

	std::vector<std::string> str_spacing;
	std::stringstream ss(stream_spacing);
	std::string item;
	while (std::getline(ss, item, '\t')) {
		str_spacing.push_back(item);
	}

	if (str_spacing.size() != 3) {
		std::cerr << "Error: spacing number should equal to 3!" << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << "Load spacing success!" << '\n' << "spacing = " << str_spacing[0] << ", " << str_spacing[1] << ", " << str_spacing[2] << std::endl;

	ImageType::SpacingType spacing;
	spacing[0] = atof(str_spacing[0].c_str());
	spacing[1] = atof(str_spacing[1].c_str());
	spacing[2] = atof(str_spacing[2].c_str());

	ImageType::SizeType resampleSize;
	for (int i = 0; i < 3; ++i) {
		resampleSize[i] = static_cast<unsigned int>(physicalSize[i] / spacing[i] + 0.5);
	}
	std::cout << "Output size: " << resampleSize << std::endl;
	std::cout << "Resampled Spacing = " << spacing[0] << ", " << spacing[1] << ", " << spacing[2] << std::endl;

	// 3. upsample
	typedef itk::IdentityTransform<double, Dimension>								TransformType;

	//typedef itk::BSplineInterpolateImageFunction<ImageType, double, double>			InterpolatorType;
	//typedef itk::LinearInterpolateImageFunction< ImageType, double >				InterpolatorType;
	typedef itk::NearestNeighborInterpolateImageFunction< ImageType, double >		InterpolatorType;
	//typedef itk::LabelImageGaussianInterpolateImageFunction<ImageType, double>		InterpolatorType;

	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	//interpolator->SetSplineOrder(3);// only for Bspline interpolate

	//interpolator->SetSigma(1.0);
	//interpolator->SetAlpha(0.3);

	typedef itk::ResampleImageFilter<ImageType, ImageType>						ResampleFilterType;

	TransformType::Pointer transform = TransformType::New();
	transform->SetIdentity();

	ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();
	resampleFilter->SetTransform(transform);
	resampleFilter->SetInterpolator(interpolator);
	resampleFilter->SetOutputOrigin(reader_candidate->GetOutput()->GetOrigin());
	resampleFilter->SetOutputSpacing(spacing);
	resampleFilter->SetSize(resampleSize);
	resampleFilter->SetOutputDirection(reader_candidate->GetOutput()->GetDirection());
	resampleFilter->SetDefaultPixelValue(0);
	resampleFilter->SetInput(reader_candidate->GetOutput());

	// 4. save the upsampled image
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(output_path.c_str());
	writer->SetInput(resampleFilter->GetOutput());
	writer->Update();

	return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////
//Goal: upsample image by itself
//Step: 1. load the upsample image candidate
//		2. compute the output size
//		3. BSpline/Linear upsample 
//		4. save the upsampled image
///////////////////////////////////////////////////////////////////////////////
bool Label_UpSample_Self(const std::string input_path, const std::string output_path, const std::string spacing_path) {

	//Step: 1. load the upsample image candidate
	typedef short												PixelType;
	const unsigned int Dimension = 3;

	typedef itk::Image<PixelType, Dimension>				ImageType;
	typedef itk::ImageFileReader<ImageType>					ReaderType;

	typedef short T_WritePixel;
	typedef itk::Image<T_WritePixel, Dimension>				WriterImageType;
	typedef itk::ImageFileWriter<WriterImageType>				WriterType;

	ReaderType::Pointer reader_candidate = ReaderType::New();

	const char* imagefile_candidate = input_path.c_str();
	if (!imagefile_candidate) {
		std::cout << "Load input file error!" << std::endl;
		return EXIT_FAILURE;
	}

	reader_candidate->SetFileName(imagefile_candidate);
	reader_candidate->Update();

	std::cout << "Read " << imagefile_candidate << " DONE!" << std::endl;

	//2. compute the output size
	ImageType::Pointer image = reader_candidate->GetOutput();
	ImageType::RegionType region = image->GetLargestPossibleRegion();
	const ImageType::SizeType size = region.GetSize();
	std::cout << "Input size: " << size << std::endl;

	const ImageType::SpacingType& sp = image->GetSpacing();

	std::cout << "Original Spacing = ";
	std::cout << sp[0] << ", " << sp[1] << ", " << sp[2] << std::endl;

	ImageType::PointType physicalSize;
	for (int i = 0; i < 3; ++i) {
		physicalSize[i] = size[i] * sp[i];
	}

	/*ImageType::SpacingType spacing;
	spacing[0] = spacing_value;
	spacing[1] = spacing_value;
	spacing[2] = spacing_value;*/
	const char* spacingfile = spacing_path.c_str();

	std::ifstream infile;
	infile.open(spacingfile, std::ifstream::in);

	if (!infile) {
		std::cerr << "Load input landmark error!" << std::endl;
		return EXIT_FAILURE;
	}

	std::string stream_spacing;
	getline(infile, stream_spacing);

	infile.close();

	std::vector<std::string> str_spacing;
	std::stringstream ss(stream_spacing);
	std::string item;
	while (std::getline(ss, item, '\t')) {
		str_spacing.push_back(item);
	}

	if (str_spacing.size() != 3) {
		std::cerr << "Error: spacing number should equal to 3!" << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << "Load spacing success!" << '\n' << "spacing = " << str_spacing[0] << ", " << str_spacing[1] << ", " << str_spacing[2] << std::endl;

	ImageType::SpacingType spacing;
	spacing[0] = atof(str_spacing[0].c_str());
	spacing[1] = atof(str_spacing[1].c_str());
	spacing[2] = atof(str_spacing[2].c_str());

	ImageType::SizeType resampleSize;
	for (int i = 0; i < 3; ++i) {
		resampleSize[i] = static_cast<unsigned int>(physicalSize[i] / spacing[i] + 0.5);
	}
	std::cout << "Output size: " << resampleSize << std::endl;
	std::cout << "Resampled Spacing = " << spacing[0] << ", " << spacing[1] << ", " << spacing[2] << std::endl;

	// 3. upsample
	typedef itk::IdentityTransform<double, Dimension>								TransformType;

	//typedef itk::BSplineInterpolateImageFunction<ImageType, double, double>			InterpolatorType;
	//typedef itk::LinearInterpolateImageFunction< ImageType, double >				InterpolatorType;
	//typedef itk::NearestNeighborInterpolateImageFunction< ImageType, double >		InterpolatorType;
	typedef itk::LabelImageGaussianInterpolateImageFunction<ImageType, double>		InterpolatorType;

	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	//interpolator->SetSplineOrder(3);// only for Bspline interpolate

	interpolator->SetSigma(1.0);
	interpolator->SetAlpha(0.3);

	typedef itk::ResampleImageFilter<ImageType, ImageType>						ResampleFilterType;

	TransformType::Pointer transform = TransformType::New();
	transform->SetIdentity();

	ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();
	resampleFilter->SetTransform(transform);
	resampleFilter->SetInterpolator(interpolator);
	resampleFilter->SetOutputOrigin(reader_candidate->GetOutput()->GetOrigin());
	resampleFilter->SetOutputSpacing(spacing);
	resampleFilter->SetSize(resampleSize);
	resampleFilter->SetOutputDirection(reader_candidate->GetOutput()->GetDirection());
	resampleFilter->SetDefaultPixelValue(0);
	resampleFilter->SetInput(reader_candidate->GetOutput());

	// 4. save the upsampled image
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(output_path.c_str());
	writer->SetInput(resampleFilter->GetOutput());
	writer->Update();

	return EXIT_SUCCESS;
}

bool Resample(const std::string input_path, const std::string output_path) {

	//Step: 1. load the upsample image candidate
	typedef short												PixelType;
	const unsigned int Dimension = 3;

	typedef itk::Image<PixelType, Dimension>				ImageType;
	typedef itk::ImageFileReader<ImageType>					ReaderType;

	//typedef short T_WritePixel;
	typedef itk::Image<PixelType, Dimension>				WriterImageType;
	typedef itk::ImageFileWriter<WriterImageType>				WriterType;

	ReaderType::Pointer reader_candidate = ReaderType::New();

	const char* imagefile_candidate = input_path.c_str();
	if (!imagefile_candidate) {
		std::cout << "Load input file error!" << std::endl;
		return EXIT_FAILURE;
	}

	reader_candidate->SetFileName(imagefile_candidate);
	reader_candidate->Update();

	std::cout << "Read " << imagefile_candidate << " DONE!" << std::endl;

	//2. compute the output size
	ImageType::Pointer image = reader_candidate->GetOutput();
	ImageType::RegionType region = image->GetLargestPossibleRegion();
	const ImageType::SizeType size = region.GetSize();
	std::cout << "Input size: " << size << std::endl;

	const ImageType::SpacingType& sp = image->GetSpacing();

	std::cout << "Original Spacing = ";
	std::cout << sp[0] << ", " << sp[1] << ", " << sp[2] << std::endl;

	ImageType::PointType physicalSize;
	for (int i = 0; i < 3; ++i) {
		physicalSize[i] = size[i] * sp[i];
	}

	ImageType::SpacingType spacing;
	spacing[0] = 2.0;
	spacing[1] = 2.0;
	spacing[2] = 2.0;
	//const char* spacingfile = spacing_path.c_str();

	//std::ifstream infile;
	//infile.open(spacingfile, std::ifstream::in);

	//if (!infile) {
	//	std::cerr << "Load input landmark error!" << std::endl;
	//	return EXIT_FAILURE;
	//}

	//std::string stream_spacing;
	//getline(infile, stream_spacing);

	//infile.close();

	//std::vector<std::string> str_spacing;
	//std::stringstream ss(stream_spacing);
	//std::string item;
	//while (std::getline(ss, item, '\t')) {
	//	str_spacing.push_back(item);
	//}

	//if (str_spacing.size() != 3) {
	//	std::cerr << "Error: spacing number should equal to 3!" << std::endl;
	//	return EXIT_FAILURE;
	//}

	//std::cout << "Load spacing success!" << '\n' << "spacing = " << str_spacing[0] << ", " << str_spacing[1] << ", " << str_spacing[2] << std::endl;

	//ImageType::SpacingType spacing;
	//spacing[0] = atof(str_spacing[0].c_str());
	//spacing[1] = atof(str_spacing[1].c_str());
	//spacing[2] = atof(str_spacing[2].c_str());

	ImageType::SizeType resampleSize;
	for (int i = 0; i < 3; ++i) {
		resampleSize[i] = static_cast<unsigned int>(physicalSize[i] / spacing[i] + 0.5);
	}
	std::cout << "Output size: " << resampleSize << std::endl;
	std::cout << "Resampled Spacing = " << spacing[0] << ", " << spacing[1] << ", " << spacing[2] << std::endl;

	resampleSize[0] = 100;
	resampleSize[1] = 100;
	resampleSize[2] = 100;

	// 3. upsample
	typedef itk::AffineTransform<double, Dimension>								TransformType;

	//typedef itk::BSplineInterpolateImageFunction<ImageType, double, TYPE>			InterpolatorType;
	//typedef itk::LinearInterpolateImageFunction< ImageType, double >				InterpolatorType;
	typedef itk::NearestNeighborInterpolateImageFunction< ImageType, double >		InterpolatorType;

	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	//interpolator->SetSplineOrder(3);// only for Bspline interpolate

	typedef itk::ResampleImageFilter<ImageType, ImageType>						ResampleFilterType;

	TransformType::Pointer transform = TransformType::New();
	transform->SetIdentity();

	//TransformType::OutputVectorType translation;
	//translation[0] = 30;
	//translation[1] = 50;
	//translation[2] = -70;
	//transform->Translate(translation);

	ImageType::PointType origin;
	origin[0] = -114;
	origin[1] = -50;
	origin[2] = 92;

	ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();
	resampleFilter->SetTransform(transform);
	resampleFilter->SetInterpolator(interpolator);
	resampleFilter->SetOutputOrigin(origin);
	resampleFilter->SetOutputSpacing(spacing);
	resampleFilter->SetSize(resampleSize);
	resampleFilter->SetOutputDirection(reader_candidate->GetOutput()->GetDirection());
	resampleFilter->SetDefaultPixelValue(50);
	resampleFilter->SetInput(reader_candidate->GetOutput());

	// 4. save the upsampled image
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(output_path.c_str());
	writer->SetInput(resampleFilter->GetOutput());
	writer->Update();

	return EXIT_SUCCESS;
}
///////////////////////////////////////////////////////////////////////////////
//Goal: resample the ct heart data 
//Step: 1. load the image 
//      2. compute the output image size
//      3. distribute the parameters for resample filter
//      4. save the resampled image
///////////////////////////////////////////////////////////////////////////////
bool DownSample(const std::string input_path, double spacing_value, const std::string output_path) {
	
	typedef short                         PixelType;
	const unsigned int                    Dimension = 3;

	typedef itk::Image<PixelType, Dimension>        ImageType;

	//Step 1. load the image from argv[1]
	typedef itk::ImageFileReader<ImageType>         ReaderType;

	ReaderType::Pointer reader = ReaderType::New();

	const char* imagefile = input_path.c_str();
	if (!imagefile) {
		std::cout << "Load input file error!" << std::endl;
		return EXIT_FAILURE;
	}

	reader->SetFileName(imagefile);
	reader->Update();

	std::cout << "Read " << imagefile << " DONE!" << std::endl;

	//Step 2. compute the output image size
	ImageType::Pointer image = reader->GetOutput();
	ImageType::RegionType region = image->GetLargestPossibleRegion();
	const ImageType::SizeType size = region.GetSize();
	std::cout << "Input size: " << size << std::endl;

	const ImageType::SpacingType& sp = image->GetSpacing();

	std::cout << "Original Spacing = ";
	std::cout << sp[0] << ", " << sp[1] << ", " << sp[2] << std::endl;

	ImageType::PointType physicalSize;
	for (int i = 0; i < 3; ++i) {
		physicalSize[i] = size[i] * sp[i];
	}

	ImageType::SpacingType spacing;
	spacing[0] = spacing_value;
	spacing[1] = spacing_value;
	spacing[2] = spacing_value;

	ImageType::SizeType resampleSize;
	for (int i = 0; i < 3; ++i) {
		resampleSize[i] = static_cast<unsigned int>(physicalSize[i] / spacing[i] + 0.5);
	}
	//std::cout << "Output size: " << resampleSize << std::endl;
	std::cout << "Resampled Spacing = " << spacing[0] << ", " << spacing[1] << ", " << spacing[2] << std::endl;

	//Step 3. distribute the parameters for resample filter

	//Beacuse this program is for downsample, you don't need interpolate like using in upsampling.
	//typedef itk::LinearInterpolateImageFunction<ImageType, TYPE>     InterpolatorType;
	typedef itk::NearestNeighborInterpolateImageFunction< ImageType, double >		InterpolatorType;
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	
	typedef itk::ResampleImageFilter<ImageType, ImageType>             ResampleFilterType;
	ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();

	typedef itk::IdentityTransform<double, 3>						     TransformType;
	TransformType::Pointer transform = TransformType::New();
	transform->SetIdentity();

	resampleFilter->SetInterpolator(interpolator);
	resampleFilter->SetSize(resampleSize);
	resampleFilter->SetOutputOrigin(image->GetOrigin());
	resampleFilter->SetOutputSpacing(spacing);
	resampleFilter->SetOutputDirection(image->GetDirection());
	resampleFilter->SetTransform(transform);
	resampleFilter->SetInput(image);
	resampleFilter->UpdateLargestPossibleRegion();
	resampleFilter->SetDefaultPixelValue(0);

	//Step 4. save the resampled image to argv[2]
	ImageType::Pointer outputImage = resampleFilter->GetOutput();
	std::cout << "Output size: " << outputImage->GetLargestPossibleRegion().GetSize() << std::endl;

	typedef itk::ImageFileWriter<ImageType> WriteType;
	WriteType::Pointer writer = WriteType::New();

	writer->SetFileName(output_path.c_str());
	writer->SetInput(outputImage);
	writer->Update();

	std::cout << "Save " << output_path << " DONE!" << std::endl;

	return EXIT_SUCCESS;
}

//GetVariance
template <typename T>
double GetVariance(std::vector<T>& data) {

	double sum = 0;

	for (auto it = data.cbegin(); it != data.cend(); ++it) {
		sum += *it;
	}
	double mean = sum / static_cast<double>(data.size());

	double varience = 0;
	for (auto it = data.cbegin(); it != data.cend(); ++it) {
		varience += (*it - mean) * (*it - mean);
	}

	return varience;
}

//GetMean
template <typename T>
double GetMean(std::vector<T>& data) {

	double sum = 0;

	for (auto it = data.cbegin(); it != data.cend(); ++it) {
		sum += *it;
	}
	double mean = sum / static_cast<double>(data.size());

	return mean;
}

///////////////////////////////////////////////////////////////////////////////
//Goal: Genegrate mask image
//Step: 1. read input image
//		2. compute the mask image //if (fg_value > threshold) {fg_iter.Set(1);
//		3. output the mask image
///////////////////////////////////////////////////////////////////////////////
void GenegrateMaskImage(const std::string input_file, const double threshold, const std::string output_file) {
	
	//1. read input image
	const unsigned int Dimension = 3;
	typedef short											PixelType;
	typedef itk::Image<PixelType, Dimension>				ImageType;
	typedef itk::ImageFileReader<ImageType>					ReaderType;

	ReaderType::Pointer reader = ReaderType::New();

	const char* imagefile = input_file.c_str();
	if (!imagefile) {
		std::cout << "Load input file error!" << std::endl;
		return;
	}

	reader->SetFileName(imagefile);
	reader->Update();

	//ImageType* image = ReadImage<PixelType>(input_file.c_str());

	std::cout << "Read " << input_file << " DONE!" << std::endl;

	//2. compute the mask image 
	ImageType::Pointer image = reader->GetOutput();

	typedef itk::ImageRegionIterator<ImageType> IteratorType;

	IteratorType fg_iter(image, image->GetRequestedRegion());

	for (fg_iter.GoToBegin(); !fg_iter.IsAtEnd(); ++fg_iter) {

		ImageType::PixelType fg_value = fg_iter.Get();

		if (fg_value > threshold) {
			fg_iter.Set(1);
		}
		else {
			fg_iter.Set(0);
		}

	}
	std::cout << "Compute the mask image DONE!" << std::endl;

	//3. output the mask image
	typedef itk::ImageFileWriter<ImageType> WriteType;
	WriteType::Pointer writer = WriteType::New();

	writer->SetFileName(output_file.c_str());
	writer->SetInput(image);
	writer->Update();

	std::cout << "Save the mask image " << output_file << " DONE!" << std::endl;

}

///////////////////////////////////////////////////////////////////////////////
//Goal: Convert landmark's the world coordinate to voxel coordinate.
//Step: 1. read the world landmark file and image data.
//      2. display the spacing, origin, scale, direction.
//      3. display the voxel coordinate and judge inside the image or not.
///////////////////////////////////////////////////////////////////////////////
bool World2Voxel(const std::string img_path, const std::string landmark_path) {
	typedef short                         PixelType;
	const unsigned int                    Dimension = 3;

	typedef itk::Image<PixelType, Dimension>      ImageType;

	// Step 1: read the world landmark file and image data.
	typedef itk::ImageFileReader<ImageType>       ReaderType;

	ReaderType::Pointer reader = ReaderType::New();

	const char* imagefile = img_path.c_str();
	if (!imagefile) {
		std::cout << "Load input file error!" << std::endl;
		return EXIT_FAILURE;
	}
	reader->SetFileName(imagefile);
	reader->Update();

	std::cout << "Read image success!" << std::endl;

	const char* landmarkfile = landmark_path.c_str();

	std::ifstream infile;
	infile.open(landmarkfile, std::ifstream::in);

	if (!infile) {
		std::cerr << "Load input landmark error!" << std::endl;
		return EXIT_FAILURE;
	}

	std::string str_landmark;
	getline(infile, str_landmark);

	infile.close();

	std::vector<std::string> str_worldLM;
	std::stringstream ss(str_landmark);
	std::string item;
	while (std::getline(ss, item, '\t')) {
		str_worldLM.push_back(item);
	}

	if (str_worldLM.size() != 3) {
		std::cerr << "Error: landmark dimension should equal to 3!" << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << "Load landmark success!" << '\n' << "Landmark's world coordinate = " << str_worldLM[0] << ", " << str_worldLM[1] << ", " << str_worldLM[2] << std::endl;

	//Step 2: display the spacing, origin, scale, direction.
	ImageType::Pointer image = reader->GetOutput();
	const ImageType::SpacingType& sp = image->GetSpacing();

	std::cout << "Spacing = ";
	std::cout << sp[0] << ", " << sp[1] << ", " << sp[2] << std::endl;

	const ImageType::PointType& origin = image->GetOrigin();

	std::cout << "Origin = ";
	std::cout << origin[0] << ", " << origin[1] << ", " << origin[2] << std::endl;

	const ImageType::RegionType scale = image->GetLargestPossibleRegion().GetSize();

	//  std::cout << "Scale = ";
	std::cout << scale << std::endl;

	const ImageType::DirectionType& direct = image->GetDirection();

	std::cout << "Direction = " << std::endl;
	std::cout << direct << std::endl;

	//Step 3: display the voxel coordinate and judge inside the image or not.
	typedef itk::Point<double, ImageType::ImageDimension> PointType;

	PointType point;

	point[0] = atof(str_worldLM[0].c_str());
	point[1] = atof(str_worldLM[1].c_str());
	point[2] = atof(str_worldLM[2].c_str());

	ImageType::IndexType pixelIndex;
	const bool isInside = image->TransformPhysicalPointToIndex(point, pixelIndex);

	std::cout << "Landmark's voxel coordinate = " << pixelIndex[0] << ", " << pixelIndex[1] << ", " << pixelIndex[2] << std::endl;

	if (isInside) {
		std::cout << "This landmark is inside the image." << std::endl;
	}
	else {
		std::cout << "This landmark is NOT inside the image." << std::endl;
	}

	return EXIT_SUCCESS;
}


///////////////////////////////////////////////////////////////////////////////
//Goal: Convert landmark's the voxel coordinate to world coordinate.
//Step: 1. read the voxel landmark file and image data.
//      2. display the spacing, origin, scale, direction.
//      3. display the world coordinate and judge inside the image or not.
///////////////////////////////////////////////////////////////////////////////
bool Voxel2World(const std::string input_path, const std::string landmark_path) {
	typedef short                         PixelType;
	const unsigned int                    Dimension = 3;

	typedef itk::Image<PixelType, Dimension>      ImageType;

	// Step 1: read the voxel landmark file and image data.
	typedef itk::ImageFileReader<ImageType>       ReaderType;

	ReaderType::Pointer reader = ReaderType::New();

	const char* imagefile = input_path.c_str();
	if (!imagefile) {
		std::cout << "Load input file error!" << std::endl;
		return EXIT_FAILURE;
	}
	reader->SetFileName(imagefile);
	reader->Update();

	std::cout << "Read image success!" << std::endl;


	const char* landmarkfile = landmark_path.c_str();

	std::ifstream infile;
	infile.open(landmarkfile, std::ifstream::in);

	if (!infile) {
		std::cerr << "Load input landmark error!" << std::endl;
		return EXIT_FAILURE;
	}

	std::string str_landmark;
	getline(infile, str_landmark);

	infile.close();

	std::vector<std::string> str_voxelLM;
	std::stringstream ss(str_landmark);
	std::string item;
	while (std::getline(ss, item, '\t')) {
		str_voxelLM.push_back(item);
	}

	if (str_voxelLM.size() != 3) {
		std::cerr << "Error: landmark dimension should equal to 3!" << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << "Load landmark success!" << '\n' << "Landmark's voxel coordinate = " << str_voxelLM[0] << ", " << str_voxelLM[1] << ", " << str_voxelLM[2] << std::endl;

	//Step 2: display the spacing, origin, scale, direction.
	ImageType::Pointer image = reader->GetOutput();
	const ImageType::SpacingType& sp = image->GetSpacing();

	std::cout << "Spacing = ";
	std::cout << sp[0] << ", " << sp[1] << ", " << sp[2] << std::endl;

	const ImageType::PointType& origin = image->GetOrigin();

	std::cout << "Origin = ";
	std::cout << origin[0] << ", " << origin[1] << ", " << origin[2] << std::endl;

	const ImageType::RegionType scale = image->GetLargestPossibleRegion().GetSize();

	//  std::cout << "Scale = ";
	std::cout << scale << std::endl;

	const ImageType::DirectionType& direct = image->GetDirection();

	std::cout << "Direction = " << std::endl;
	std::cout << direct << std::endl;

	//Step 3: display the world coordinate and judge inside the image or not.
	typedef itk::Point<double, ImageType::ImageDimension> PointType;

	ImageType::IndexType point;

	point[0] = atoi(str_voxelLM[0].c_str());
	point[1] = atoi(str_voxelLM[1].c_str());
	point[2] = atoi(str_voxelLM[2].c_str());

	PointType worldLM;
	image->TransformIndexToPhysicalPoint(point, worldLM);

	std::cout << "Landmark's world coordinate = " << worldLM[0] << ", " << worldLM[1] << ", " << worldLM[2] << std::endl;

	int isInside = -1;
	if (point[0] < 0 || point[1] < 0 || point[2] < 0 ||
		point[0] > image->GetLargestPossibleRegion().GetSize()[0] ||
		point[1] > image->GetLargestPossibleRegion().GetSize()[1] ||
		point[2] > image->GetLargestPossibleRegion().GetSize()[2]) {
		isInside = 0;
	}
	else {
		isInside = 1;
	}

	if (isInside == 1) {
		std::cout << "This landmark is inside the image." << std::endl;
	}
	else if (isInside == 0) {
		std::cout << "This landmark is NOT inside the image." << std::endl;
	}
	else {
		std::cout << "isInside error!" << std::endl;
	}

	return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////
//Goal: Histogram matching for input image based on reference image
//Step: 1. load the input image and reference image
//      2. config the histogram matching
//      3. output the matched image to file
///////////////////////////////////////////////////////////////////////////////
bool HistogramMatching(const std::string ref_path, const std::string input_path, const std::string output_path) {
	
	//1. load the input image and reference image
	typedef short                         PixelType;
	const unsigned int                    Dimension = 3;

	typedef itk::Image<PixelType, Dimension>      ImageType;

	typedef itk::ImageFileReader<ImageType>       ReaderType;

	//load input image
	ReaderType::Pointer reader = ReaderType::New();

	const char* imagefile = input_path.c_str();
	if (!imagefile) {
		std::cout << "Load input file error!" << std::endl;
		return EXIT_FAILURE;
	}
	reader->SetFileName(imagefile);
	reader->Update();

	std::cout << "Read " << imagefile << " DONE!" << std::endl;

	//load reference image 
	ReaderType::Pointer reader_refer = ReaderType::New();

	const char* imagefile_refer = ref_path.c_str();
	if (!imagefile_refer) {
		std::cout << "Load input file error!" << std::endl;
		return EXIT_FAILURE;
	}

	reader_refer->SetFileName(imagefile_refer);
	reader_refer->Update();

	std::cout << "Read " << imagefile_refer << " DONE!" << std::endl;

	//2. config the histogram matching
	typedef float                                      InternalPixelType;
	typedef itk::Image< InternalPixelType, Dimension > InternalImageType;
	typedef itk::CastImageFilter< ImageType, InternalImageType >		InputImageCasterType;
	typedef itk::CastImageFilter< ImageType, InternalImageType >		ReferImageCasterType;
	
	InputImageCasterType::Pointer inputImageCaster = InputImageCasterType::New();
	ReferImageCasterType::Pointer referImageCaster = ReferImageCasterType::New();
	
	inputImageCaster->SetInput(reader->GetOutput());
	referImageCaster->SetInput(reader_refer->GetOutput());

	typedef itk::HistogramMatchingImageFilter<InternalImageType, InternalImageType >		 MatchingFilterType;
	MatchingFilterType::Pointer matcher = MatchingFilterType::New();

	matcher->SetInput(inputImageCaster->GetOutput());
	matcher->SetReferenceImage(referImageCaster->GetOutput());

	matcher->SetNumberOfHistogramLevels(1024);
	matcher->SetNumberOfMatchPoints(7);

	matcher->ThresholdAtMeanIntensityOn();

	//3. output the matched image to file
	typedef itk::CastImageFilter<InternalImageType, ImageType >          CastFilterType;
	typedef itk::ImageFileWriter< ImageType >  WriterType;
	
	WriterType::Pointer      writer = WriterType::New();
	CastFilterType::Pointer  caster = CastFilterType::New();
	
	writer->SetFileName(output_path.c_str());
	caster->SetInput(matcher->GetOutput());
	writer->SetInput(caster->GetOutput());
	writer->Update();

	std::cout << "Write " << output_path << " DONE!" << std::endl;

	return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////
//Goal: generate the segmentation based on different label's probability map
//Step: 1. load foreground probability map
//		2. compare the value of probability map with threshold and set the label to label map
//      3. output the label map to file
///////////////////////////////////////////////////////////////////////////////
bool GenegrateSegmentation(const std::string input_path, const double threshold, const std::string output_path) {
	typedef double            		  PixelType;

	const unsigned int                  Dimension = 3;

	typedef itk::Image<PixelType, Dimension>      ImageType;

	//1. load two probability map
	typedef itk::ImageFileReader<ImageType>       ReaderType;

	//foreground_image
	ReaderType::Pointer reader = ReaderType::New();

	const char* filename = input_path.c_str();
	if (!filename) {
		std::cout << "Load input file error!" << std::endl;
		return EXIT_FAILURE;
	}
	reader->SetFileName(filename);
	reader->Update();

	std::cout << "Read foreground probability map success!" << std::endl;

	//2. compare the value of probability map with threshold and set the label to label map
	ImageType::Pointer foreground_image = reader->GetOutput();

	typedef itk::ImageRegionIterator<ImageType> IteratorType;

	IteratorType fg_iter(foreground_image, foreground_image->GetRequestedRegion());

	for (fg_iter.GoToBegin(); !fg_iter.IsAtEnd(); ++fg_iter) {

		ImageType::PixelType fg_value = fg_iter.Get();

		if (fg_value >= threshold) {
			fg_iter.Set(1);
		}
		else {
			fg_iter.Set(0);
		}

	}

	//Step 3: output the label map to file
	//typedef itk::Image<PixelType, Dimension>      OutputImageType;
	typedef itk::ImageFileWriter<ImageType> WriteType;
	WriteType::Pointer writer = WriteType::New();

	writer->SetFileName(output_path.c_str());
	writer->SetInput(foreground_image);
	writer->Update();

	std::cout << "Save the segmenation " << output_path << " DONE!" << std::endl;

	return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////
//Goal: Extract largest connected component from label image
//Step: 1. load the input image
//      2. compute the largest connected component
//      3. output the largest connected component to file
///////////////////////////////////////////////////////////////////////////////
bool ExtractLargestComponent(const std::string input_path, const std::string output_path) {

	const unsigned int Dimension = 3;
	typedef unsigned short					  PixelType;

	//1. load the input image
	typedef itk::Image< PixelType, Dimension >  ImageType;
	typedef itk::ImageFileReader< ImageType >   ReaderType;

	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(input_path.c_str());

	typedef itk::Image< PixelType, Dimension > OutputImageType;
	typedef itk::ConnectedComponentImageFilter <ImageType, OutputImageType >   ConnectedComponentImageFilterType;

	ConnectedComponentImageFilterType::Pointer connected = ConnectedComponentImageFilterType::New();
	connected->SetInput(reader->GetOutput());
	connected->Update();

	std::cout << "Number of objects: " << connected->GetObjectCount() << std::endl;

	//2. compute the largest connected component
	typedef itk::LabelShapeKeepNObjectsImageFilter< OutputImageType > LabelShapeKeepNObjectsImageFilterType;
	LabelShapeKeepNObjectsImageFilterType::Pointer labelShapeKeepNObjectsImageFilter = LabelShapeKeepNObjectsImageFilterType::New();
	labelShapeKeepNObjectsImageFilter->SetInput(connected->GetOutput());
	labelShapeKeepNObjectsImageFilter->SetBackgroundValue(0);
	labelShapeKeepNObjectsImageFilter->SetNumberOfObjects(1);
	labelShapeKeepNObjectsImageFilter->SetAttribute(LabelShapeKeepNObjectsImageFilterType::LabelObjectType::NUMBER_OF_PIXELS);

	typedef itk::RescaleIntensityImageFilter< OutputImageType, ImageType > RescaleFilterType;
	RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
	rescaleFilter->SetOutputMinimum(0);
	rescaleFilter->SetOutputMaximum(1);
	rescaleFilter->SetInput(labelShapeKeepNObjectsImageFilter->GetOutput());

	//3. output the largest connected component to file
	ImageType::Pointer outputImage = rescaleFilter->GetOutput();

	typedef itk::ImageFileWriter<OutputImageType> WriteType;
	WriteType::Pointer writer = WriteType::New();

	writer->SetFileName(output_path.c_str());
	writer->SetInput(outputImage);
	writer->Update();

	std::cout << "Save " << output_path << " DONE!" << std::endl;

	return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////
//Goal: fill the holes in the label image
//Step: 1. load the label image
//      2. fill the hole by the filter
//      3. output the revised label image to file
///////////////////////////////////////////////////////////////////////////////
bool HillHole(const std::string input_path, const std::string output_path, int radius) {
	//1. load the label image
	typedef   unsigned short         InputPixelType;
	typedef   unsigned short         OutputPixelType;
	typedef itk::Image< InputPixelType, 3 >   InputImageType;
	typedef itk::Image< OutputPixelType, 3 >   OutputImageType;

	typedef itk::ImageFileReader< InputImageType  >  ReaderType;
	typedef itk::ImageFileWriter< OutputImageType >  WriterType;
	ReaderType::Pointer reader = ReaderType::New();
	WriterType::Pointer writer = WriterType::New();
	reader->SetFileName(input_path.c_str());
	writer->SetFileName(output_path.c_str());

	InputImageType::Pointer image = reader->GetOutput();
	std::cout << "Load image success!" << std::endl;

	//2. fill the hole by the filter
	typedef itk::VotingBinaryHoleFillingImageFilter< InputImageType, OutputImageType >  FilterType;
	FilterType::Pointer filter = FilterType::New();

	const unsigned int radiusX = radius;
	const unsigned int radiusY = radiusX;
	const unsigned int radiusZ = radiusX;

	InputImageType::SizeType indexRadius;
	indexRadius[0] = radiusX; // radius along x
	indexRadius[1] = radiusY; // radius along y
	indexRadius[2] = radiusZ;
	filter->SetRadius(indexRadius);

	filter->SetBackgroundValue(0);
	filter->SetForegroundValue(1);
	filter->SetMajorityThreshold(2);
	filter->SetInput(image);
	writer->SetInput(filter->GetOutput());
	std::cout << "Fill hole success!" << std::endl;
	// 3. output the revised label image to file
	writer->Update();

	std::cout << "Output to file success!" << std::endl;
	return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////
//Goal: correct the image to origin(0,0,0), direction(1,1,1), isotropy spacing
//Step: 1. load the inputted image
//      2. correct the image
//      3. output the revised image to file
///////////////////////////////////////////////////////////////////////////////
bool Correct(const std::string referenceImage, const std::string inputImageFile, const std::string outputImageFile, double spacing_value) {

	//Step: 1. load the reference image
	typedef short												PixelType;
	const unsigned int Dimension = 3;

	typedef itk::Image<PixelType, Dimension>				ImageType;
	typedef itk::ImageFileReader<ImageType>					ReaderType;

	typedef short T_WritePixel;
	typedef itk::Image<T_WritePixel, Dimension>				WriterImageType;
	typedef itk::ImageFileWriter<WriterImageType>				WriterType;

	ReaderType::Pointer reader_refer = ReaderType::New();

	const char* imagefile_refer = referenceImage.c_str();
	if (!imagefile_refer) {
		std::cout << "Load input file error!" << std::endl;
		return EXIT_FAILURE;
	}

	reader_refer->SetFileName(imagefile_refer);
	reader_refer->Update();

	std::cout << "Read " << imagefile_refer << " DONE!" << std::endl;

	ImageType::Pointer refer_image = reader_refer->GetOutput();

	//2. load the upsample image candidate
	ReaderType::Pointer reader_candidate = ReaderType::New();

	const char* imagefile_candidate = inputImageFile.c_str();
	if (!imagefile_candidate) {
		std::cout << "Load input file error!" << std::endl;
		return EXIT_FAILURE;
	}

	reader_candidate->SetFileName(imagefile_candidate);
	reader_candidate->Update();

	std::cout << "Read " << imagefile_candidate << " DONE!" << std::endl;

	ImageType::Pointer image = reader_candidate->GetOutput();

	ImageType::SpacingType spacing;
	spacing[0] = spacing_value;
	spacing[1] = spacing_value;
	spacing[2] = spacing_value;

	// 3. upsample
	typedef itk::IdentityTransform<double, Dimension>								TransformType;

	//typedef itk::BSplineInterpolateImageFunction<ImageType, double, TYPE>			InterpolatorType;
	//typedef itk::LinearInterpolateImageFunction< ImageType, double >				InterpolatorType;
	typedef itk::NearestNeighborInterpolateImageFunction< ImageType, double >		InterpolatorType;

	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	//interpolator->SetSplineOrder(3);// only for Bspline interpolate

	typedef itk::ResampleImageFilter<ImageType, ImageType>						ResampleFilterType;

	TransformType::Pointer transform = TransformType::New();
	transform->SetIdentity();

	ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();
	resampleFilter->SetTransform(transform);
	resampleFilter->SetInterpolator(interpolator);
	resampleFilter->SetOutputOrigin(image->GetOrigin());
	resampleFilter->SetOutputSpacing(image->GetSpacing());
	resampleFilter->SetSize(image->GetLargestPossibleRegion().GetSize());
	resampleFilter->SetOutputDirection(refer_image->GetDirection());
	resampleFilter->SetDefaultPixelValue(0);
	resampleFilter->SetInput(image);

	// 4. save the upsampled image
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(outputImageFile.c_str());
	writer->SetInput(resampleFilter->GetOutput());
	writer->Update();

	return EXIT_SUCCESS;

	////1. load the inputted image
	//typedef short											PixelType;
	//const unsigned int Dimension = 3;

	//typedef itk::Image<PixelType, Dimension>				ImageType;
	//typedef itk::ImageFileReader<ImageType>					ReaderType;

	//typedef short											WritePixelType;
	//typedef itk::Image<WritePixelType, Dimension>			WriterImageType;
	//typedef itk::ImageFileWriter<WriterImageType>			WriterType;

	//ReaderType::Pointer reader = ReaderType::New();

	//const char* imagefile = input_path.c_str();
	//if (!imagefile) {
	//	std::cout << "Load input file error!" << std::endl;
	//	return EXIT_FAILURE;
	//}

	//reader->SetFileName(imagefile);
	//reader->Update();

	//std::cout << "Read " << imagefile << " DONE!" << std::endl;
	//
	////2. correct the image

	//ImageType::Pointer image = reader->GetOutput();
	//ImageType::RegionType region = image->GetLargestPossibleRegion();
	//const ImageType::SizeType size = region.GetSize();
	//std::cout << "Original size = " << size << std::endl;

	//const ImageType::SpacingType& sp = image->GetSpacing();

	//std::cout << "Original Spacing = ";
	//std::cout << sp[0] << ", " << sp[1] << ", " << sp[2] << std::endl;

	//ImageType::PointType physicalSize;
	//for (int i = 0; i < 3; ++i) {
	//	physicalSize[i] = size[i] * sp[i];
	//}

	//const ImageType::PointType& origin = image->GetOrigin();

	//std::cout << "Original Origin = ";
	//std::cout << origin[0] << ", " << origin[1] << ", " << origin[2] << std::endl;

	//const ImageType::DirectionType& direct = image->GetDirection();

	//std::cout << "Original Direction = " << std::endl;
	//std::cout << direct << std::endl;

	////revise spacing
	//ImageType::SpacingType newSpacing;
	//newSpacing[0] = spacing_value;
	//newSpacing[1] = spacing_value;
	//newSpacing[2] = spacing_value;

	////revise size
	//ImageType::SizeType newSize;
	//for (int i = 0; i < 3; ++i) {
	//	newSize[i] = static_cast<unsigned int>(physicalSize[i] / newSpacing[i] + 0.5);
	//}

	//std::cout << "New Size = " << newSize << std::endl;
	//std::cout << "New Spacing = " << newSpacing[0] << ", " << newSpacing[1] << ", " << newSpacing[2] << std::endl;

	////revise origin
	//ImageType::PointType newOrigin;
	//newOrigin.Fill(0.0);

	//std::cout << "New Origin = ";
	//std::cout << newOrigin[0] << ", " << newOrigin[1] << ", " << newOrigin[2] << std::endl;

	////revise direction
	//ImageType::DirectionType newDirection;
	//newDirection.SetIdentity();

	//std::cout << "New Direction = " << std::endl;
	//std::cout << newDirection << std::endl;

	//ImageType::Pointer inputImage = reader->GetOutput();

	//typedef itk::ChangeInformationImageFilter< ImageType >			FilterType;
	//FilterType::Pointer filter = FilterType::New();

	//filter->SetOutputSpacing(newSpacing);
	//filter->ChangeSpacingOn();

	//filter->SetOutputOrigin(newOrigin);
	//filter->ChangeOriginOn();

	//filter->SetOutputDirection(newDirection);
	//filter->ChangeDirectionOn();
	//
	////3. output the revised image to file
	//WriterType::Pointer writer = WriterType::New();
	//writer->SetFileName(output_path.c_str());
	//writer->SetInput(filter->GetOutput());
	//writer->Update();

	//return EXIT_SUCCESS;

	/*
	//Step: 1. load the inputted image
	typedef short											PixelType;
	const unsigned int Dimension = 3;

	typedef itk::Image<PixelType, Dimension>				ImageType;
	typedef itk::ImageFileReader<ImageType>					ReaderType;

	typedef short											WritePixelType;
	typedef itk::Image<WritePixelType, Dimension>			WriterImageType;
	typedef itk::ImageFileWriter<WriterImageType>			WriterType;

	ReaderType::Pointer reader = ReaderType::New();

	const char* imagefile = input_path.c_str();
	if (!imagefile) {
		std::cout << "Load input file error!" << std::endl;
		return EXIT_FAILURE;
	}

	reader->SetFileName(imagefile);
	reader->Update();

	std::cout << "Read " << imagefile << " DONE!" << std::endl;


	//2. correct the image
	
	ImageType::Pointer image = reader->GetOutput();
	ImageType::RegionType region = image->GetLargestPossibleRegion();
	const ImageType::SizeType size = region.GetSize();
	std::cout << "Original size = " << size << std::endl;

	const ImageType::SpacingType& sp = image->GetSpacing();

	std::cout << "Original Spacing = ";
	std::cout << sp[0] << ", " << sp[1] << ", " << sp[2] << std::endl;

	ImageType::PointType physicalSize;
	for (int i = 0; i < 3; ++i) {
		physicalSize[i] = size[i] * sp[i];
	}

	const ImageType::PointType& origin = image->GetOrigin();

	std::cout << "Original Origin = ";
	std::cout << origin[0] << ", " << origin[1] << ", " << origin[2] << std::endl;

	const ImageType::DirectionType& direct = image->GetDirection();

	std::cout << "Original Direction = " << std::endl;
	std::cout << direct << std::endl;

	//revise spacing
	ImageType::SpacingType newSpacing;
	newSpacing[0] = spacing_value;
	newSpacing[1] = spacing_value;
	newSpacing[2] = spacing_value;

	//revise size
	ImageType::SizeType newSize;
	for (int i = 0; i < 3; ++i) {
		newSize[i] = static_cast<unsigned int>(physicalSize[i] / newSpacing[i] + 0.5);
	}

	std::cout << "New Size = " << newSize << std::endl;
	std::cout << "New Spacing = " << newSpacing[0] << ", " << newSpacing[1] << ", " << newSpacing[2] << std::endl;

	//revise origin
	ImageType::PointType newOrigin;
	newOrigin.Fill(0.0);
	
	std::cout << "New Origin = ";
	std::cout << newOrigin[0] << ", " << newOrigin[1] << ", " << newOrigin[2] << std::endl;

	//revise direction
	ImageType::DirectionType newDirection;
	newDirection.SetIdentity();

	std::cout << "New Direction = " << std::endl;
	std::cout << newDirection << std::endl;

	// set linear interpolate method
	typedef itk::IdentityTransform<double, Dimension>								TransformType;

	//typedef itk::BSplineInterpolateImageFunction<ImageType, double, TYPE>       InterpolatorType;

	typedef itk::LinearInterpolateImageFunction< ImageType, double >				InterpolatorType;
	InterpolatorType::Pointer interpolator = InterpolatorType::New();

	typedef itk::ResampleImageFilter<ImageType, ImageType>						ResampleFilterType;

	TransformType::Pointer transform = TransformType::New();
	transform->SetIdentity();

	//InterpolatorType::Pointer interpolator = InterpolatorType::New();
	//interpolator->SetSplineOrder(3);

	ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();
	resampleFilter->SetInput(reader->GetOutput());
	resampleFilter->SetTransform(transform);
	//resampleFilter->SetInterpolator(interpolator);
	resampleFilter->SetOutputOrigin(newOrigin);
	resampleFilter->SetOutputSpacing(newSpacing);
	resampleFilter->SetSize(newSize);
	resampleFilter->SetOutputDirection(newDirection);
	resampleFilter->UpdateLargestPossibleRegion();
	resampleFilter->SetDefaultPixelValue(0);
	

	// 3. output the revised image to file
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(output_path.c_str());
	writer->SetInput(resampleFilter->GetOutput());
	writer->Update();

	//WriteImage<PixelType, short>(resampleFilter->GetOutput(), outputImageFile.c_str());


	return EXIT_SUCCESS;

	*/
}

///////////////////////////////////////////////////////////////////////////////
//Goal: dilate the binary image
//Step: 1. load the binary image
//      2. config the dilate filter
//      3. output the dilated binary image to file
///////////////////////////////////////////////////////////////////////////////
bool Dilate(const std::string input_path, const std::string output_path, unsigned int radius) {
	
	const int Dimension = 3;
	typedef unsigned char					PixelType;
	//1. load the binary image
	typedef itk::Image<PixelType, Dimension>    ImageType;
	typedef itk::ImageFileReader<ImageType> ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(input_path.c_str());

	//2. config the dilate filter
	PixelType foreground = 1;
	typedef itk::BinaryBallStructuringElement<ImageType::PixelType, 3>      StructuringElementType;
	StructuringElementType structuringElement;
	structuringElement.SetRadius(radius); // radius = 1 ==> 3x3 structuring element
	structuringElement.CreateStructuringElement();

	typedef itk::BinaryDilateImageFilter <ImageType, ImageType, StructuringElementType> 	BinaryDilateImageFilterType;

	BinaryDilateImageFilterType::Pointer dilateFilter = BinaryDilateImageFilterType::New();
	dilateFilter->SetDilateValue(foreground);
	dilateFilter->SetInput(reader->GetOutput());
	dilateFilter->SetKernel(structuringElement);

	//3. output the dilated binary image to file
	typedef unsigned char									WritePixelType;
	typedef itk::Image<WritePixelType, Dimension>			WriterImageType;
	typedef itk::ImageFileWriter<WriterImageType>			WriterType;

	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(output_path.c_str());
	writer->SetInput(dilateFilter->GetOutput());
	writer->Update();

	return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////
//Goal: dilate the erode image
//Step: 1. load the erode image
//      2. config the erode filter
//      3. output the eroded binary image to file
///////////////////////////////////////////////////////////////////////////////
bool Erode(const std::string input_path, const std::string output_path, unsigned int radius) {

	const int Dimension = 3;
	typedef unsigned char					PixelType;
	//1. load the binary image
	typedef itk::Image<PixelType, Dimension>    ImageType;
	typedef itk::ImageFileReader<ImageType> ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(input_path.c_str());

	//2. config the erode filter
	PixelType foreground = 1;
	typedef itk::BinaryBallStructuringElement<ImageType::PixelType, 3>      StructuringElementType;
	StructuringElementType structuringElement;
	structuringElement.SetRadius(radius); // radius = 1 ==> 3x3 structuring element
	structuringElement.CreateStructuringElement();

	typedef itk::BinaryErodeImageFilter <ImageType, ImageType, StructuringElementType> 	BinaryErodeImageFilterType;

	BinaryErodeImageFilterType::Pointer erodeFilter = BinaryErodeImageFilterType::New();
	erodeFilter->SetErodeValue(foreground);
	erodeFilter->SetInput(reader->GetOutput());
	erodeFilter->SetKernel(structuringElement);

	//3. output the eroded binary image to file
	typedef unsigned char									WritePixelType;
	typedef itk::Image<WritePixelType, Dimension>			WriterImageType;
	typedef itk::ImageFileWriter<WriterImageType>			WriterType;

	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(output_path.c_str());
	writer->SetInput(erodeFilter->GetOutput());
	writer->Update();

	return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////
//Goal: convert data type to 'nii.gz'
//Step: 1. load the input image
//      2. output the image to nii.gz file
///////////////////////////////////////////////////////////////////////////////
bool Convert(const std::string input_path, const std::string output_path) {
	
	const int Dimension = 3;
	typedef short								PixelType;
	//1. load the input image
	typedef itk::Image<PixelType, Dimension>    ImageType;
	typedef itk::ImageFileReader<ImageType>		ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(input_path.c_str());

	//2. output the image to nii.gz file
	typedef short											WritePixelType;
	typedef itk::Image<WritePixelType, Dimension>			WriterImageType;
	typedef itk::ImageFileWriter<WriterImageType>			WriterType;

	typedef itk::CastImageFilter< ImageType, WriterImageType > CastFilterType;
	CastFilterType::Pointer castFilter = CastFilterType::New();
	castFilter->SetInput( reader->GetOutput() );

	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(output_path.c_str());
	writer->SetInput( castFilter->GetOutput() );
	writer->Update();

	return EXIT_SUCCESS;

}

///////////////////////////////////////////////////////////////////////////////
//Goal: count voxel number
//Step: 1. load the input image
//      2. count the voxel number using point iterator
///////////////////////////////////////////////////////////////////////////////
bool Count(const std::string input_path) {
	
	typedef unsigned short              PixelType;
	const unsigned int                  Dimension = 3;

	typedef itk::Image<PixelType, Dimension>      ImageType;


	//1. load the input image
	typedef itk::ImageFileReader<ImageType>       ReaderType;

	ReaderType::Pointer reader = ReaderType::New();


	reader->SetFileName(input_path.c_str());
	reader->Update();

	std::cout << "Read image success!" << std::endl;

	//2. count the voxel number using point iterator
	ImageType::Pointer image = reader->GetOutput();

	unsigned int all_num = 0, foreg_num = 0, backg_num = 0;

	//  typedef itk::ImageRegionConstIterator<ImageType> ConstIteratorType;
	typedef itk::ImageRegionIterator<ImageType> IteratorType;

	IteratorType iter(image, image->GetRequestedRegion());

	for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {

		ImageType::PixelType value = iter.Get();

		if (value == 0) {
			backg_num++;
		}
		else if (value == 1) {
			foreg_num++;
		}
		else {
			std::cerr << "Label is " << value << "?" << std::endl;
		}

		all_num++;
	}

	std::cout << "Foreground voxel number = " << foreg_num << "/" << all_num << " = " << (double)foreg_num / all_num << std::endl;
	std::cout << "Background voxel number = " << backg_num << "/" << all_num << " = " << (double)backg_num / all_num << std::endl;

	return EXIT_SUCCESS;

}

///////////////////////////////////////////////////////////////////////////////
//Goal: get the bounding box from binary image
//Step: 1. read the binary image
//		2. get the bounding box
//		3. transform the start point and end point from voxel coordinate to world coordinate
//		4. output the two points world coordinate to txt file
///////////////////////////////////////////////////////////////////////////////
bool BoundingBox(const std::string input_path, const std::string start_pt, const std::string end_pt) {

	//1. read the binary image
	typedef itk::ImageMaskSpatialObject<3>      ImageMaskSpatialObject;
	typedef ImageMaskSpatialObject::ImageType   ImageType;
	typedef ImageType::RegionType               RegionType;
	typedef itk::ImageFileReader< ImageType >   ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(input_path.c_str());
	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject & excp)
	{
		std::cerr << excp << std::endl;
		return EXIT_FAILURE;
	}

	//2. get the bounding box
	ImageType::Pointer image = reader->GetOutput();

	ImageMaskSpatialObject::Pointer maskSO = ImageMaskSpatialObject::New();
	maskSO->SetImage(image);
	RegionType boundingBoxRegion = maskSO->GetAxisAlignedBoundingBoxRegion();
	std::cout << "Bounding Box Region: " << boundingBoxRegion << std::endl;

	ImageType::IndexType box_start = boundingBoxRegion.GetIndex();
	ImageType::SizeType box_size = boundingBoxRegion.GetSize();

	std::cout << "Box start index = " << box_start[0] << " " << box_start[1] << " " << box_start[2] << "\n";
	std::cout << "Box size = " << box_size[0] << " " << box_size[1] << " " << box_size[2] << "\n";

	ImageType::IndexType box_end;
	box_end[0] = box_start[0] + box_size[0] - 1;
	box_end[1] = box_start[1] + box_size[1] - 1;
	box_end[2] = box_start[2] + box_size[2] - 1;

	std::cout << "Box end index = " << box_end[0] << " " << box_end[1] << " " << box_end[2] << "\n";

	//3. transform the start point and end point from voxel coordinate to world coordinate
	typedef itk::Point<double, ImageType::ImageDimension> PointType;
	PointType world_start, world_end;

	image->TransformIndexToPhysicalPoint(box_start, world_start);
	image->TransformIndexToPhysicalPoint(box_end, world_end);

	std::cout << "start point's world coordinate = " << world_start[0] << ", " << world_start[1] << ", " << world_start[2] << std::endl;
	std::cout << "end point's world coordinate = " << world_end[0] << ", " << world_end[1] << ", " << world_end[2] << std::endl;

	//4. output the two points world coordinate to txt file
	//start point 
	std::ofstream start_file(start_pt, std::ios::out);
	if (!start_file)
	{
		std::cout << "Error opening file for writing." << std::endl;
		return EXIT_FAILURE;
	}

	start_file << world_start[0] << " " << world_start[1] << " " << world_start[2] << "\n";

	start_file.close();

	//end point
	std::ofstream end_file(end_pt, std::ios::out);
	if (!end_file)
	{
		std::cout << "Error opening file for writing." << std::endl;
		return EXIT_FAILURE;
	}

	end_file << world_end[0] << " " << world_end[1] << " " << world_end[2] << "\n";

	end_file.close();

	return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////
//Goal: crop volume
//Step: 1. load the input image, start and end point world coordinate
//      2. crop the volume
//		3. output the image to file
///////////////////////////////////////////////////////////////////////////////
bool Crop(const std::string input_path, const std::string start, const std::string end, const int pad, const std::string output_path) {
	
	typedef short              PixelType;
	const unsigned int                  Dimension = 3;

	typedef itk::Image<PixelType, Dimension>      ImageType;

	//1. load the input image, start and end point world coordinate
	typedef itk::ImageFileReader<ImageType>       ReaderType;

	ReaderType::Pointer reader = ReaderType::New();

	reader->SetFileName(input_path.c_str());
	reader->Update();

	std::cout << "Read image success!" << std::endl;

	ImageType::Pointer image = reader->GetOutput();

	//start point world coordinate
	const char* startfile = start.c_str();

	std::ifstream infile_start;
	infile_start.open(startfile, std::ifstream::in);

	if (!infile_start) {
		std::cerr << "Load input startfile error!" << std::endl;
		return EXIT_FAILURE;
	}

	std::string str_startindex;
	getline(infile_start, str_startindex);

	infile_start.close();

	std::vector<std::string> vector_start;
	std::stringstream ss(str_startindex);
	std::string item;
	while (std::getline(ss, item, ' ')) {
		vector_start.push_back(item);
	}

	if (vector_start.size() != 3) {
		std::cerr << "Error: vector_start dimension should equal to 3!" << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << "Load start world coordinate success!" << '\n' << "start point's world coordinate = "
		<< vector_start[0] << ", " << vector_start[1] << ", " << vector_start[2] << std::endl;

	//end point world coordinate
	const char* endfile = end.c_str();

	std::ifstream infile_end;
	infile_end.open(endfile, std::ifstream::in);

	if (!infile_end) {
		std::cerr << "Load input endfile error!" << std::endl;
		return EXIT_FAILURE;
	}

	std::string str_end;
	getline(infile_end, str_end);

	infile_end.close();

	std::vector<std::string> vector_end;
	std::stringstream ss2(str_end);
	std::string item2;
	while (std::getline(ss2, item2, ' ')) {
		vector_end.push_back(item2);
	}

	if (vector_end.size() != 3) {
		std::cerr << "Error: vector_end dimension should equal to 3!" << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << "Load vector_end success!" << '\n' << "end point's world coordinate = " 
		<< vector_end[0] << ", " << vector_end[1] << ", " << vector_end[2] << std::endl;

	//2. crop the volume
	std::cout << "Image largest region: " << image->GetLargestPossibleRegion() << std::endl;
	ImageType::SizeType largestSize = image->GetLargestPossibleRegion().GetSize();

	//transform the start end point from world coordinate to voxel coordinate
	typedef itk::Point<double, ImageType::ImageDimension> PointType;

	PointType start_point;

	start_point[0] = atof(vector_start[0].c_str());
	start_point[1] = atof(vector_start[1].c_str());
	start_point[2] = atof(vector_start[2].c_str());

	PointType end_point;

	end_point[0] = atof(vector_end[0].c_str());
	end_point[1] = atof(vector_end[1].c_str());
	end_point[2] = atof(vector_end[2].c_str());

	ImageType::IndexType startIndex, endIndex;
	const bool startInside = image->TransformPhysicalPointToIndex(start_point, startIndex);
	const bool endInside = image->TransformPhysicalPointToIndex(end_point, endIndex);

	std::cout << "start point's index = " << startIndex[0] << ", " << startIndex[1] << ", " << startIndex[2] << std::endl;
	std::cout << "end point's index = " << endIndex[0] << ", " << endIndex[1] << ", " << endIndex[2] << std::endl;

	if (startInside) {
		std::cout << "This startInside is inside the image." << std::endl;
	}
	else {
		std::cout << "This startInside is NOT inside the image." << std::endl;
	}

	if (endInside) {
		std::cout << "This endInside is inside the image." << std::endl;
	}
	else {
		std::cout << "This endInside is NOT inside the image." << std::endl;
	}

	//convert physical size to image size
	ImageType::SpacingType spacing = image->GetSpacing();
	std::cout << "spacing = " << spacing[0] << ", " << spacing[1] << ", " << spacing[2] << std::endl;

	int padmm[3];
	for (int i = 0; i < 3; i++) {
		padmm[i] = static_cast<int>(pad / spacing[i] + 0.5);
	}

	std::cout << "pad mm = " << padmm[0] << ", " << padmm[1] << ", " << padmm[2] << std::endl;

	//revise the bounding box size
	for (int i = 0; i < 3; i++) {
		int temp = startIndex[i] - padmm[i];
		if (temp < 0) {
			startIndex[i] = 0;
		}
		else
			startIndex[i] = temp;
	}

	for (int i = 0; i < 3; i++) {
		int temp2 = endIndex[i] + padmm[i];
		if (temp2 > largestSize[i]) {
			endIndex[i] = largestSize[i];
		}
		else
			endIndex[i] = temp2;
	}

	ImageType::SizeType desiredSize;
	for (int i = 0; i < 3; i++) {
		desiredSize[i] = endIndex[i] - startIndex[i];
	}
	
	ImageType::RegionType desiredRegion(startIndex, desiredSize);

	std::cout << "desiredRegion: " << desiredRegion << std::endl;

	typedef itk::ExtractImageFilter< ImageType, ImageType > FilterType;
	FilterType::Pointer filter = FilterType::New();
	filter->SetExtractionRegion(desiredRegion);
	filter->SetInput(image);
#if ITK_VERSION_MAJOR >= 4
	filter->SetDirectionCollapseToIdentity(); // This is required.
#endif
	filter->Update();

	//3. output the image to file
	typedef itk::Image<PixelType, Dimension>				WriterImageType;
	typedef itk::ImageFileWriter<WriterImageType>			WriterType;

	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(output_path.c_str());
	writer->SetInput(filter->GetOutput());
	writer->Update();

	return EXIT_SUCCESS;

}

///////////////////////////////////////////////////////////////////////////////
//Goal: crop volume for femur
//Step: 1. load the input label image, lesser trochanter landmark point world coordinate
//      2. tansfer the world coordinate to image index
//		3. crop the volume based on lesser trochanter image index
//		4. output the image to file
///////////////////////////////////////////////////////////////////////////////
bool CropFemur(const std::string input_path, const std::string landmark_path, const std::string output_path) {
	
	// 1. load the input label image, lesser trochanter landmark point world coordinate
	typedef short                         PixelType;
	const unsigned int                    Dimension = 3;

	typedef itk::Image<PixelType, Dimension>      ImageType;

	typedef itk::ImageFileReader<ImageType>       ReaderType;

	ReaderType::Pointer reader = ReaderType::New();

	const char* imagefile = input_path.c_str();
	if (!imagefile) {
		std::cout << "Load input file error!" << std::endl;
		return EXIT_FAILURE;
	}
	reader->SetFileName(imagefile);
	reader->Update();

	std::cout << "Read image success!" << std::endl;

	//lesser trochanter landmark point world coordinate
	const char* landmarkfile = landmark_path.c_str();

	std::ifstream infile;
	infile.open(landmarkfile, std::ifstream::in);

	if (!infile) {
		std::cerr << "Load input landmark error!" << std::endl;
		return EXIT_FAILURE;
	}

	std::string str_landmark;
	getline(infile, str_landmark);

	infile.close();

	std::vector<std::string> str_worldLM;
	std::stringstream ss(str_landmark);
	std::string item;
	while (std::getline(ss, item, '\t')) {
		str_worldLM.push_back(item);
	}

	if (str_worldLM.size() != 3) {
		std::cerr << "Error: landmark dimension should equal to 3!" << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << "Load landmark success!" << '\n' << "Landmark's world coordinate = " << str_worldLM[0] << ", " << str_worldLM[1] << ", " << str_worldLM[2] << std::endl;

	// //Step 2: display the spacing, origin, scale, direction.
	 ImageType::Pointer image = reader->GetOutput();
	// const ImageType::SpacingType& sp = image->GetSpacing();

	// std::cout << "Spacing = ";
	// std::cout << sp[0] << ", " << sp[1] << ", " << sp[2] << std::endl;

	// const ImageType::PointType& origin = image->GetOrigin();

	// std::cout << "Origin = ";
	// std::cout << origin[0] << ", " << origin[1] << ", " << origin[2] << std::endl;

	const ImageType::SizeType scale = image->GetLargestPossibleRegion().GetSize();

	//  std::cout << "Scale = ";
	std::cout << scale << std::endl;

	// const ImageType::DirectionType& direct = image->GetDirection();

	// std::cout << "Direction = " << std::endl;
	// std::cout << direct << std::endl;

	//2. tansfer the world coordinate to image index
	typedef itk::Point<double, ImageType::ImageDimension> PointType;

	PointType point;

	point[0] = atof(str_worldLM[0].c_str());
	point[1] = atof(str_worldLM[1].c_str());
	point[2] = atof(str_worldLM[2].c_str());

	ImageType::IndexType pixelIndex;
	const bool isInside = image->TransformPhysicalPointToIndex(point, pixelIndex);

	std::cout << "Landmark's voxel coordinate = " << pixelIndex[0] << ", " << pixelIndex[1] << ", " << pixelIndex[2] << std::endl;

	if (isInside) {
		std::cout << "This landmark is inside the image." << std::endl;
	}
	else {
		std::cout << "This landmark is NOT inside the image." << std::endl;
	}

	//3. crop the volume based on lesser trochanter image index
	ImageType::SizeType regionSize;
	regionSize[0] = scale[0];
	regionSize[1] = scale[1];
	regionSize[2] = pixelIndex[2];

	
	ImageType::IndexType  regionIndex;
	for (int i=0; i<3; i++) {
		regionIndex[i] = 0;
	}

	ImageType::RegionType region;
	region.SetSize(regionSize);
	region.SetIndex(regionIndex);

	itk::ImageRegionIterator<ImageType> imageIterator(image, region);

	while(!imageIterator.IsAtEnd())
	{
	// Get the value of the current pixel
	//unsigned char val = imageIterator.Get();
	//std::cout << (int)val << std::endl;

	// Set the current pixel to white
	imageIterator.Set(0);

	++imageIterator;
	}


	//output the image to file
	typedef itk::Image<PixelType, Dimension>				WriterImageType;
	typedef itk::ImageFileWriter<WriterImageType>			WriterType;

	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(output_path.c_str());
	writer->SetInput(image);
	writer->Update();

	return EXIT_SUCCESS;

}

///////////////////////////////////////////////////////////////////////////////
//Goal: image minus a specified value
//Step: 1. load the input image
//      2. minus the specified value iteratorly
//		3. output the image to file
///////////////////////////////////////////////////////////////////////////////
bool Minus(const std::string input_path, const int minus_value, const std::string output_path) {
	typedef short					    PixelType;
	const unsigned int                  Dimension = 3;

	typedef itk::Image<PixelType, Dimension>      ImageType;

	//1. load the input image, start and end point world coordinate
	typedef itk::ImageFileReader<ImageType>       ReaderType;

	ReaderType::Pointer reader = ReaderType::New();

	reader->SetFileName(input_path.c_str());
	reader->Update();

	std::cout << "Read image success!" << std::endl;

	ImageType::Pointer image = reader->GetOutput();

	//2. minus the specified value iteratorly
	typedef itk::ImageRegionIterator<ImageType> IteratorType;

	IteratorType iter(image, image->GetRequestedRegion());

	for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {

		ImageType::PixelType value = iter.Get();

		iter.Set(value - minus_value);
	
	}

	//3. output the image to file
	typedef itk::Image<PixelType, Dimension>				WriterImageType;
	typedef itk::ImageFileWriter<WriterImageType>			WriterType;

	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(output_path.c_str());
	writer->SetInput(image);
	writer->Update();

	return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////
//Goal: put cropped image back to original image
//Step: 1. load the reference and input image
//      2. create and initialize the new image same size as reference image
//		3. transform the world coordinate to the voxel coordinate
//		4. assign the segmentation to new image by image iterator
//		5. output the image to file
///////////////////////////////////////////////////////////////////////////////
bool PutBack(const std::string ref_path, const std::string input_path, const std::string start_pt, const std::string end_pt, const std::string output_path) {
	
	//----------   1. load the reference and input image   ----------  
	typedef short					    PixelType;
	const unsigned int                  Dimension = 3;

	typedef itk::Image<PixelType, Dimension>      ImageType;

	typedef itk::ImageFileReader<ImageType>       ReaderType;

	//ref image
	ReaderType::Pointer ref_reader = ReaderType::New();

	ref_reader->SetFileName(ref_path.c_str());
	ref_reader->Update();

	std::cout << "Read ref image success!" << std::endl;

	ImageType::Pointer ref_image = ref_reader->GetOutput();

	//input image
	ReaderType::Pointer in_reader = ReaderType::New();

	in_reader->SetFileName(input_path.c_str());
	in_reader->Update();

	std::cout << "Read input image success!" << std::endl;

	ImageType::Pointer in_image = in_reader->GetOutput();

	//----------   2. create and initialize the new image same size as reference image   ----------  
	ImageType::Pointer new_image = ImageType::New();
	new_image->SetOrigin(ref_image->GetOrigin());
	new_image->SetSpacing(ref_image->GetSpacing());
	new_image->SetRegions(ref_image->GetLargestPossibleRegion());
	new_image->SetDirection(ref_image->GetDirection());
	new_image->Allocate();
	new_image->FillBuffer(0);//initialize 0 for each voxel

	//----------   3. transform the world coordinate to the voxel coordinate   ----------  
	//load start point world coordinate
	const char* startfile = start_pt.c_str();

	std::ifstream infile_start;
	infile_start.open(startfile, std::ifstream::in);

	if (!infile_start) {
		std::cerr << "Load input startfile error!" << std::endl;
		return EXIT_FAILURE;
	}

	std::string str_startindex;
	getline(infile_start, str_startindex);

	infile_start.close();

	std::vector<std::string> vector_start;
	std::stringstream ss(str_startindex);
	std::string item;
	while (std::getline(ss, item, ' ')) {
		vector_start.push_back(item);
	}

	if (vector_start.size() != 3) {
		std::cerr << "Error: vector_start dimension should equal to 3!" << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << "\nLoad start world coordinate success!" << '\n' << "start point's world coordinate = "
		<< vector_start[0] << ", " << vector_start[1] << ", " << vector_start[2] << std::endl;

	//load end point world coordinate
	const char* endfile = end_pt.c_str();

	std::ifstream infile_end;
	infile_end.open(endfile, std::ifstream::in);

	if (!infile_end) {
		std::cerr << "Load input endfile error!" << std::endl;
		return EXIT_FAILURE;
	}

	std::string str_end;
	getline(infile_end, str_end);

	infile_end.close();

	std::vector<std::string> vector_end;
	std::stringstream ss2(str_end);
	std::string item2;
	while (std::getline(ss2, item2, ' ')) {
		vector_end.push_back(item2);
	}

	if (vector_end.size() != 3) {
		std::cerr << "Error: vector_end dimension should equal to 3!" << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << "\nLoad vector_end success!" << '\n' << "end point's world coordinate = "
		<< vector_end[0] << ", " << vector_end[1] << ", " << vector_end[2] << std::endl;


	//transform world to voxel coordinate
	typedef itk::Point<double, ImageType::ImageDimension> PointType;

	PointType start_point;

	start_point[0] = atof(vector_start[0].c_str());
	start_point[1] = atof(vector_start[1].c_str());
	start_point[2] = atof(vector_start[2].c_str());

	PointType end_point;

	end_point[0] = atof(vector_end[0].c_str());
	end_point[1] = atof(vector_end[1].c_str());
	end_point[2] = atof(vector_end[2].c_str());

	ImageType::IndexType startIndex, endIndex;
	const bool startInside = ref_image->TransformPhysicalPointToIndex(start_point, startIndex);
	const bool endInside = ref_image->TransformPhysicalPointToIndex(end_point, endIndex);

	std::cout << "start point's index = " << startIndex[0] << ", " << startIndex[1] << ", " << startIndex[2] << std::endl;
	std::cout << "end point's index = " << endIndex[0] << ", " << endIndex[1] << ", " << endIndex[2] << std::endl;

	if (startInside) {
		std::cout << "This startInside is inside the image." << std::endl;
	}
	else {
		std::cout << "This startInside is NOT inside the image." << std::endl;
		return EXIT_FAILURE;
	}

	if (endInside) {
		std::cout << "This endInside is inside the image." << std::endl;
	}
	else {
		std::cout << "This endInside is NOT inside the image." << std::endl;
		return EXIT_FAILURE;
	}

	//revise the bounding box size
	for (int i = 0; i < 3; i++) {
		int temp = startIndex[i] - 40;
		if (temp < 0) {
			startIndex[i] = 0;
		}
		else
			startIndex[i] = temp;
	}

	ImageType::SizeType largestSize = ref_image->GetLargestPossibleRegion().GetSize();

	for (int i = 0; i < 3; i++) {
		int temp2 = endIndex[i] + 40;
		if (temp2 > largestSize[i]) {
			endIndex[i] = largestSize[i];
		}
		else
			endIndex[i] = temp2;
	}

	//----------   4. assign the segmentation to new image by image iterator   ----------  
	ImageType::SizeType in_size = in_image->GetLargestPossibleRegion().GetSize();
	std::cout << "input size = " << in_size << std::endl;

	ImageType::SizeType crop_size;
	for (int i = 0; i < 3; i++) {
		crop_size[i] = endIndex[i] - startIndex[i] + 1;
	}
	std::cout << "crop size = " << crop_size << std::endl;

	for (int i = 0; i < 3; i++) {
		if (in_size[i] != crop_size[i]) {
			std::cerr << "Error: Size is not euqal!" << std::endl;
			return EXIT_FAILURE;
		}
	}

	typedef itk::ImageRegionConstIterator<ImageType>	ConstIteratorType;
	typedef itk::ImageRegionIterator<ImageType>			IteratorType;

	ImageType::RegionType out_region;
	out_region.SetIndex(startIndex);
	out_region.SetSize(in_size);

	IteratorType new_iter(new_image, out_region);
	ConstIteratorType in_iter(in_image, in_image->GetRequestedRegion());

	new_iter.GoToBegin();
	in_iter.GoToBegin();

	while ( !in_iter.IsAtEnd() ) {

		new_iter.Set( in_iter.Get() );

		++new_iter;
		++in_iter;
	}


	//----------   5. output the image to file   ----------  
	typedef itk::Image<PixelType, Dimension>				WriterImageType;
	typedef itk::ImageFileWriter<WriterImageType>			WriterType;

	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName(output_path.c_str());
	writer->SetInput(new_image);
	writer->Update();

	return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////
//Goal: label transform
//Step: 1. Read the CT data image
//      2. Change the old label value to 0-7.
//		3. Save the CT data image to NiFTI type data 
///////////////////////////////////////////////////////////////////////////////
bool LabelTransform(const std::string input_path, const std::string output_path) {
	typedef unsigned short						PixelType;
	const unsigned int                  Dimension = 3;

	typedef itk::Image<PixelType, Dimension>      ImageType;

	//Step 1: Read the CT data image
	typedef itk::ImageFileReader<ImageType>       ReaderType;

	ReaderType::Pointer reader = ReaderType::New();

	const char* filename = input_path.c_str();
	if (!filename) {
		std::cout << "Load input file error!" << std::endl;
		return EXIT_FAILURE;
	}
	reader->SetFileName(filename);
	reader->Update();

	std::cout << "Read image success!" << std::endl;

	//Step 2: Change the old label value to 0-7.
	ImageType::Pointer image = reader->GetOutput();

	//  typedef itk::ImageRegionConstIterator<ImageType> ConstIteratorType;
	typedef itk::ImageRegionIterator<ImageType> IteratorType;

	IteratorType iter(image, image->GetRequestedRegion());

	int count = 0;
	for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {

		ImageType::PixelType value = iter.Get();

		if (value == 0) {
			continue;
		}
		else if (value == 500) {
			iter.Set(1);
		}
		else if (value == 600) {
			iter.Set(2);
		}
		else if (value == 420 || value == 421) {
			iter.Set(3);
		}
		else if (value == 550) {
			iter.Set(4);
		}
		else if (value == 205) {
			iter.Set(5);
		}
		else if (value == 820) {
			iter.Set(6);
		}
		else if (value == 850) {
			iter.Set(7);
		}
		else {
			std::cerr << "Error: other label value !" << "\t Value = " << value << std::endl;
		}

		//count++;
	}
	//std::cout << "count = " << count << std::endl;
	std::cout << "Change the label success!" << std::endl;

	//Step 3: Save the CT data image to NiFTI type data 
	typedef itk::ImageFileWriter<ImageType> WriteType;
	WriteType::Pointer writer = WriteType::New();

	writer->SetFileName(output_path.c_str());
	writer->SetInput(image);
	writer->Update();

	std::cout << "Write transformed image to file success!" << std::endl;

	return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////
//Goal: adaptive histogram equalization
//Step: 1. load the input image
//      2. adaptive histogram equalization
//		3. output the image to file
///////////////////////////////////////////////////////////////////////////////
bool AdaptiveHistogramEqualization(const std::string input_path, const std::string output_path, const double alpha, const double beta, unsigned int radius) {

	//1. load the input image
	typedef short						PixelType;
	const unsigned int                  Dimension = 3;

	typedef itk::Image<PixelType, Dimension>      ImageType;

	typedef itk::ImageFileReader<ImageType>       ReaderType;

	ReaderType::Pointer reader = ReaderType::New();

	const char* filename = input_path.c_str();
	if (!filename) {
		std::cout << "Load input file error!" << std::endl;
		return EXIT_FAILURE;
	}
	reader->SetFileName(filename);
	reader->Update();

	std::cout << "Read image success!" << std::endl;

	ImageType::Pointer image = reader->GetOutput();

	// 2. adaptive histogram equalization
	typedef  itk::AdaptiveHistogramEqualizationImageFilter< ImageType > AdaptiveHistogramEqualizationImageFilterType;
	AdaptiveHistogramEqualizationImageFilterType::Pointer filter = 
		AdaptiveHistogramEqualizationImageFilterType::New();
	
	filter->SetInput(image);
	filter->SetAlpha(alpha);
	filter->SetBeta(beta);
	filter->SetRadius(radius);

	filter->Update();

	//3. output the image to file
	typedef itk::ImageFileWriter<ImageType> WriteType;
	WriteType::Pointer writer = WriteType::New();

	writer->SetFileName(output_path.c_str());
	writer->SetInput(filter->GetOutput());
	writer->Update();

	std::cout << "Write transformed image to file success!" << std::endl;

	return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////
//Goal: rescale intensity
//Step: 1. load the input image
//      2. rescale image intensity to [min, max]
//		3. output the image to file
///////////////////////////////////////////////////////////////////////////////
bool RescaleIntensity(const std::string input_path, unsigned int min, unsigned int max, const std::string output_path) {

	//1. load the input image
	typedef short						PixelType;
	const unsigned int                  Dimension = 3;

	typedef itk::Image<PixelType, Dimension>      ImageType;

	typedef itk::ImageFileReader<ImageType>       ReaderType;

	ReaderType::Pointer reader = ReaderType::New();

	const char* filename = input_path.c_str();
	if (!filename) {
		std::cout << "Load input file error!" << std::endl;
		return EXIT_FAILURE;
	}
	reader->SetFileName(filename);
	reader->Update();

	std::cout << "Read image success!" << std::endl;

	ImageType::Pointer image = reader->GetOutput();

	//2. rescale image intensity to [min, max]
	typedef itk::RescaleIntensityImageFilter< ImageType, ImageType > RescaleFilterType;
	RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
	rescaleFilter->SetInput(image);
	rescaleFilter->SetOutputMinimum(min);
	rescaleFilter->SetOutputMaximum(max);

	//3. output the image to file
	typedef itk::ImageFileWriter<ImageType> WriteType;
	WriteType::Pointer writer = WriteType::New();

	writer->SetFileName(output_path.c_str());
	writer->SetInput(rescaleFilter->GetOutput());
	writer->Update();

	std::cout << "Write transformed image to file success!" << std::endl;

	return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////
//Goal: orient image to specified direction
//Step: 1. load the input image
//      2. orient image to specified direction
//		3. output the image to file
///////////////////////////////////////////////////////////////////////////////
bool Orient(const std::string input_path, const std::string output_path) {
	
	//1. load the input image
	typedef short						PixelType;
	const unsigned int                  Dimension = 3;

	typedef itk::Image<PixelType, Dimension>      ImageType;

	typedef itk::ImageFileReader<ImageType>       ReaderType;

	ReaderType::Pointer reader = ReaderType::New();

	//To use the deprecated classes (AnalyzeImageIO) you must turn ITKV3_COMPATIBILITY:BOOL=ON during configuration. 
	//itk::AnalyzeImageIO::Pointer io = itk::AnalyzeImageIO::New();

	const char* filename = input_path.c_str();
	if (!filename) {
		std::cout << "Load input file error!" << std::endl;
		return EXIT_FAILURE;
	}

	//reader->SetImageIO(io);
	reader->SetFileName(filename);
	reader->Update();

	std::cout << "Read image success!" << std::endl;

	ImageType::Pointer image = reader->GetOutput();

	//2. orient image to specified direction
	itk::OrientImageFilter<ImageType, ImageType>::Pointer orienter = itk::OrientImageFilter<ImageType, ImageType>::New();
	orienter->UseImageDirectionOff();
	orienter->SetGivenCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_AIL);
	orienter->SetDesiredCoordinateOrientation(itk::SpatialOrientation::ITK_COORDINATE_ORIENTATION_AIR);
	orienter->SetInput(image);
	orienter->Update();
	
	//3. output the image to file
	typedef itk::ImageFileWriter<ImageType> WriteType;
	WriteType::Pointer writer = WriteType::New();

	writer->SetFileName(output_path.c_str());
	writer->SetInput(orienter->GetOutput());
	writer->Update();

	std::cout << "Write transformed image to file success!" << std::endl;

	return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////
//Goal: judge two image is same or not
//Step: 1. load the reference image and input image 
//      2. compare every voxel of two image and print the first 1000 value.
//		3. output the conclusion
///////////////////////////////////////////////////////////////////////////////
bool Judge(const std::string referenceImage, const std::string input_path) {
	
	//1. load the reference image and input image
	typedef short											PixelType;
	const unsigned int Dimension = 3;

	typedef itk::Image<PixelType, Dimension>				ImageType;
	typedef itk::ImageFileReader<ImageType>					ReaderType;

	ReaderType::Pointer reader_refer = ReaderType::New();

	const char* imagefile_refer = referenceImage.c_str();
	if (!imagefile_refer) {
		std::cout << "Load input file error!" << std::endl;
		return EXIT_FAILURE;
	}

	reader_refer->SetFileName(imagefile_refer);
	reader_refer->Update();

	std::cout << "Read " << imagefile_refer << " DONE!" << std::endl;

	ImageType::Pointer ref_img = reader_refer->GetOutput();

	//input image
	ReaderType::Pointer reader_candidate = ReaderType::New();

	const char* imagefile_candidate = input_path.c_str();
	if (!imagefile_candidate) {
		std::cout << "Load input file error!" << std::endl;
		return EXIT_FAILURE;
	}

	reader_candidate->SetFileName(imagefile_candidate);
	reader_candidate->Update();

	ImageType::Pointer input_img = reader_candidate->GetOutput();

	std::cout << "Read " << imagefile_candidate << " DONE!" << std::endl;

	//2. compare every voxel of two image and print the first 1000 value.
	typedef itk::ImageRegionIterator<ImageType> IteratorType;

	IteratorType iter_ref(ref_img, ref_img->GetRequestedRegion());
	IteratorType iter(input_img, input_img->GetRequestedRegion());
	
	int count = 0;
	for (iter.GoToBegin(), iter_ref.GoToBegin(); !iter.IsAtEnd(); ++iter, ++iter_ref) {

		ImageType::PixelType value_ref = iter_ref.Get();
		ImageType::PixelType value = iter.Get();

		if (count < 1000) {
			std::cout << "value_ref = " << value_ref << "  vs.  " << "value = " << value << std::endl;
			count++;
			continue;
		}

		if (value == value_ref) {
			count++;
			continue;
		}
		else
		{
			// 3. output the conclusion
			std::cout << "Two image is different!" << std::endl;
			return EXIT_FAILURE;
		}

		count++;
	}
	// 3. output the conclusion
	std::cout << "Two image is same!" << std::endl;

	return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////
//Goal: cut off the histogram
//Step: 1. load input image 
//      2. sort all the voxel by intensity value from max to min
//		3. get the max*ratio's value and set all the larger value to the max*ratio's value
//		4. output the image to file
///////////////////////////////////////////////////////////////////////////////
bool CutOff(const std::string input_path, const std::string output_path, const double ratio) {
	
	//1. load input image
	typedef short						PixelType;
	const unsigned int                  Dimension = 3;

	typedef itk::Image<PixelType, Dimension>      ImageType;

	typedef itk::ImageFileReader<ImageType>       ReaderType;

	ReaderType::Pointer reader = ReaderType::New();

	const char* filename = input_path.c_str();
	if (!filename) {
		std::cout << "Load input file error!" << std::endl;
		return EXIT_FAILURE;
	}

	//reader->SetImageIO(io);
	reader->SetFileName(filename);
	reader->Update();

	std::cout << "Read image success!" << std::endl;

	ImageType::Pointer image = reader->GetOutput();

	//2. sort all the voxel by intensity value from max to min
	unsigned int all_num = 0, foreg_num = 0, backg_num = 0;

	typedef itk::ImageRegionConstIterator<ImageType> ConstIteratorType;
	//typedef itk::ImageRegionIterator<ImageType> IteratorType;

	ConstIteratorType iter(image, image->GetRequestedRegion());

	std::vector<short> v_intensity;

	int count = 0;
	for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {

		ImageType::PixelType value = iter.Get();

		v_intensity.push_back(value);

		count++;
	}

	if (count != v_intensity.size()) {
		std::cout << "Error: count != v_intensity.size()" << std::endl;
	}
	else {
		std::cout << "The number of all voxel is " << v_intensity.size() << std::endl;
	}
	
	std::sort(v_intensity.begin(), v_intensity.end());

	//std::cout << "v_intensity contains:";
	//for (std::vector<short>::iterator it = v_intensity.end(); it != v_intensity.begin(); it--)
	//	std::cout << ' ' << *it;
	//std::cout << '\n';

	//3. get the max*ratio's value and set all the larger value to the max*ratio's value
	short max_value = v_intensity[v_intensity.size() - 1];
	std::cout << "max intensity = " << max_value << std::endl;

	int index = static_cast<int>(v_intensity.size() * (1 - ratio));
	std::cout << "cut off index = " << index << std::endl;

	short cutoff_value = v_intensity[index - 1];
	std::cout << "cut off intensity = " << cutoff_value << std::endl;

	typedef itk::ImageRegionIterator<ImageType> IteratorType;

	IteratorType iter2(image, image->GetRequestedRegion());

	for (iter2.GoToBegin(); !iter2.IsAtEnd(); ++iter2) {

		ImageType::PixelType value = iter2.Get();

		if (value > cutoff_value) {
			iter2.Set(cutoff_value);
		}
	}

	//3. output the image to file
	typedef itk::ImageFileWriter<ImageType> WriteType;
	WriteType::Pointer writer = WriteType::New();

	writer->SetFileName(output_path.c_str());
	writer->SetInput(image);
	writer->Update();

	std::cout << "Write transformed image to file success!" << std::endl;

	return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////////////////////////
//Goal: extract left-top and right-bottom landmarks of femur
//Step: 1. read the binary image
//      2. get the bounding box and find the lt and rb index
//		3. pad the lt and rb points
//		4. output the new mask image to file
//		5. output the lt and rb world coordinates to file
///////////////////////////////////////////////////////////////////////////////
bool ExtractLandmarkOfFemur(const std::string input_path, const std::string start_pt, const std::string end_pt, const std::string output_path) {

	//1. read the binary image
	typedef itk::ImageMaskSpatialObject<3>      ImageMaskSpatialObject;
	typedef ImageMaskSpatialObject::ImageType   ImageType;
	typedef ImageType::RegionType               RegionType;
	typedef itk::ImageFileReader< ImageType >   ReaderType;
	ReaderType::Pointer reader = ReaderType::New();
	reader->SetFileName(input_path.c_str());
	try
	{
		reader->Update();
	}
	catch (itk::ExceptionObject & excp)
	{
		std::cerr << excp << std::endl;
		return EXIT_FAILURE;
	}

	//2. get the bounding box and find the lt and rb index
	ImageType::Pointer image = reader->GetOutput();

	ImageMaskSpatialObject::Pointer maskSO = ImageMaskSpatialObject::New();
	maskSO->SetImage(image);
	RegionType boundingBoxRegion = maskSO->GetAxisAlignedBoundingBoxRegion();
	std::cout << "Bounding Box Region: " << boundingBoxRegion << std::endl;

	ImageType::IndexType box_start = boundingBoxRegion.GetIndex();
	ImageType::SizeType box_size = boundingBoxRegion.GetSize();

	std::cout << "Box start index = " << box_start[0] << " " << box_start[1] << " " << box_start[2] << "\n";
	
	ImageType::IndexType box_end;
	box_end[0] = box_start[0] + box_size[0] - 1;
	box_end[1] = box_start[1] + box_size[1] - 1;
	box_end[2] = box_start[2] + box_size[2] - 1;

	std::cout << "Box end index = " << box_end[0] << " " << box_end[1] << " " << box_end[2] << "\n";
	std::cout << "Box size = " << box_size[0] << " " << box_size[1] << " " << box_size[2] << "\n";

	//3. pad the lt and rb points
	//NOTE: Need not check the lt and rb points whether inside the image volume as landmarks can outside the image volume.
	ImageType::IndexType lt, rb;
	for (int i=0; i<3; i++) {
		lt[i] = box_end[i] + static_cast<int>(box_size[i]/10.0 + 1.0);
	}

	for (int i=0; i<3; i++) {
		rb[i] = box_start[i] - static_cast<int>(box_size[i]/10.0 + 1.0);
	}

	std::cout << "left-top index = " << lt[0] << " " << lt[1] << " " << lt[2] << "\n";
	std::cout << "right-bottom index  = " << rb[0] << " " << rb[1] << " " << rb[2] << "\n";


	//4. output the new mask image to file 
	RegionType newRegion;
	ImageType::SizeType newSize;
	for (int i=0; i<3; i++) {
		newSize[i] = lt[i] - rb[i] + 1;
	}
	std::cout << "New box size = " << newSize[0] << " " << newSize[1] << " " << newSize[2] << "\n";

	newRegion.SetIndex(rb);
	newRegion.SetSize(newSize);

	itk::ImageRegionIterator<ImageType> imageIterator(image, newRegion);

	while(!imageIterator.IsAtEnd())
	{
	// Get the value of the current pixel
	//unsigned char val = imageIterator.Get();
	//std::cout << (int)val << std::endl;

	// Set the current pixel to white
	imageIterator.Set(1);

	++imageIterator;
	}

	typedef itk::ImageFileWriter<ImageType> WriteType;
	WriteType::Pointer writer = WriteType::New();

	writer->SetFileName(output_path.c_str());
	writer->SetInput(image);
	writer->Update();

	std::cout << "Write transformed image to file success!" << std::endl;

	//5. output the lt and rb world coordinates to file
	typedef itk::Point<double, ImageType::ImageDimension> PointType;
	PointType world_start, world_end;

	image->TransformIndexToPhysicalPoint(rb, world_start);
	image->TransformIndexToPhysicalPoint(lt, world_end);

	std::cout << "start point's world coordinate = " << world_start[0] << ", " << world_start[1] << ", " << world_start[2] << std::endl;
	std::cout << "end point's world coordinate = " << world_end[0] << ", " << world_end[1] << ", " << world_end[2] << std::endl;


	//start point 
	std::ofstream start_file(start_pt, std::ios::out);
	if (!start_file)
	{
		std::cout << "Error opening file for writing." << std::endl;
		return EXIT_FAILURE;
	}

	start_file << world_start[0] << "\t" << world_start[1] << "\t" << world_start[2] << "\n";

	start_file.close();

	//end point
	std::ofstream end_file(end_pt, std::ios::out);
	if (!end_file)
	{
		std::cout << "Error opening file for writing." << std::endl;
		return EXIT_FAILURE;
	}

	end_file << world_end[0] << "\t" << world_end[1] << "\t" << world_end[2] << "\n";

	end_file.close();

	return EXIT_SUCCESS;

}



int main(int argc, char * argv[])
{
	try {

		using namespace boost::program_options;

		options_description generic("Generic options");
		generic.add_options()("help", "display helper information");

		options_description functional("Functional options (must pick one)");
		functional.add_options()("nearest_upsample", "upsample an image");
		functional.add_options()("linear_upsample", "upsample an image"); 
		functional.add_options()("upsample_self", "upsample an image by itself");
		functional.add_options()("label_upsample_self", "upsample an image by label interpolation");
		functional.add_options()("resample", "resample image"); 
		functional.add_options()("downsample", "downsample an image");
		functional.add_options()("mask", "convert image format/change data type");
		functional.add_options()("world2voxel", "convert point from world coord to voxel coord");
		functional.add_options()("voxel2world", "convert point from voxel coord to world coord");
		functional.add_options()("histogram_match", "match the input image's histogram based on reference image");
		functional.add_options()("generate_segment", "generate the segmentation based on probability map");
		functional.add_options()("largest_component", "extract largest connected component from label image");
		functional.add_options()("fill_hole", "extract largest connected component from label image");
		functional.add_options()("correct", "correct the image to origin(0,0,0), direction(1,1,1), isotropy spacing");
		functional.add_options()("dilate", "dilate binary image");
		functional.add_options()("erode", "erode binary image"); 
		functional.add_options()("convert", "convert data type to 'nii.gz'"); 
		functional.add_options()("count", "count voxel number");
		functional.add_options()("crop", "crop volumn");
		functional.add_options()("crop_femur", "crop volumn for femur");
		functional.add_options()("boundingbox", "get bounding box");
		functional.add_options()("minus", "image minus specified value");
		functional.add_options()("putback", "put cropped image back to original image");	
		functional.add_options()("label_transform", "image label transform");
		functional.add_options()("adaptive_histogram_equalization", "adaptive histogram equalization");
		functional.add_options()("rescale_intensity", "rescale intensity");
		functional.add_options()("orient", "orient image to specified direction");	
		functional.add_options()("judge", "judge two image is same or not");
		functional.add_options()("cutoff", "cut off the histogram");
		functional.add_options()("extract_landmark_femur", "extract left-top and right-bottom landmarks of femur");
		
		
		options_description params("Parameter options");
		params.add_options()("in", value<std::string>(), "input image");
		params.add_options()("ref", value<std::string>(), "reference image");
		params.add_options()("out", value<std::string>(), "reference image");
		params.add_options()("threshold", value<double>(), "mask threshold");
		params.add_options()("landmark", value<std::string>(), "landmark position");
		params.add_options()("spacing", value<double>(), "spacing value");
		params.add_options()("radius", value<int>(), "neighbour radius");
		params.add_options()("pad", value<int>(), "pad");
		params.add_options()("start_pt", value<std::string>(), "start point file");
		params.add_options()("end_pt", value<std::string>(), "end point file");
		params.add_options()("minus_value", value<int>(), "minus value");
		params.add_options()("min", value<int>(), "minimum value");
		params.add_options()("max", value<int>(), "maximum value");
		params.add_options()("alpha", value<double>(), "alpha");
		params.add_options()("beta", value<double>(), "beta");
		params.add_options()("ratio", value<double>(), "ratio");
		params.add_options()("spacing_file", value<std::string>(), "spacing file");

		options_description cmd_options;
		cmd_options.add(generic).add(functional).add(params);

		variables_map vm;
		store(command_line_parser(argc, argv).options(cmd_options).style(command_line_style::unix_style ^ command_line_style::allow_short).run(), vm);
		notify(vm);

		std::cout << "Read Input Parameters DONE" << std::endl;

		if (vm.count("help") || argc == 1) //help
		{
			std::cerr << cmd_options << std::endl;
			// std::cout << "Hello, ImageTools!" << std::endl;
			print_help();
			return -1;
		}
		
		
		if (vm.count("nearest_upsample")) // nearest upsample
		{
			std::string input_path = Get<std::string>(vm, "in");
			std::string ref_path = Get<std::string>(vm, "ref");
			std::string output_path = Get<std::string>(vm, "out");
			double spacing_value = Get<double>(vm, "spacing");

			std::cout << "input_path = " << input_path << std::endl;
			std::cout << "ref_path = " << ref_path << std::endl;
			std::cout << "output_path = " << output_path << std::endl;
			std::cout << "spacing_value = " << spacing_value << std::endl;

			//Nearest_UpSample(ref_path, input_path, output_path, spacing_value);

			LoadPixelType type = ReadPixelType(input_path.c_str());

			switch (type)
			{
			case _UCHAR: Nearest_UpSample<unsigned char>(ref_path, input_path, output_path, spacing_value);  break;
			case _CHAR: Nearest_UpSample<char>(ref_path, input_path, output_path, spacing_value);  break;
			case _USHORT: Nearest_UpSample<unsigned short>(ref_path, input_path, output_path, spacing_value);  break;
			case _SHORT: Nearest_UpSample<short>(ref_path, input_path, output_path, spacing_value);  break;
			case _UINT: Nearest_UpSample<unsigned int>(ref_path, input_path, output_path, spacing_value);  break;
			case _INT: Nearest_UpSample<int>(ref_path, input_path, output_path, spacing_value);  break;
			case _ULONG: Nearest_UpSample<unsigned long>(ref_path, input_path, output_path, spacing_value);  break;
			case _LONG: Nearest_UpSample<long>(ref_path, input_path, output_path, spacing_value);  break;
			case _FLOAT: Nearest_UpSample<float>(ref_path, input_path, output_path, spacing_value);  break;
			case _DOUBLE: Nearest_UpSample<double>(ref_path, input_path, output_path, spacing_value);  break;
			default: std::cout << "unrecognized pixel type"; break;
			}
		}
		else if (vm.count("linear_upsample")) // linear upsample
		{
			std::string input_path = Get<std::string>(vm, "in");
			std::string ref_path = Get<std::string>(vm, "ref");
			std::string output_path = Get<std::string>(vm, "out");
			double spacing_value = Get<double>(vm, "spacing");

			std::cout << "input_path = " << input_path << std::endl;
			std::cout << "ref_path = " << ref_path << std::endl;
			std::cout << "output_path = " << output_path << std::endl;
			std::cout << "spacing_value = " << spacing_value << std::endl;

			//Linear_UpSample(ref_path, input_path, output_path, spacing_value);
		
			LoadPixelType type = ReadPixelType(input_path.c_str());

			switch (type)
			{
			case _UCHAR: Linear_UpSample<unsigned char>(ref_path, input_path, output_path, spacing_value);  break;
			case _CHAR: Linear_UpSample<char>(ref_path, input_path, output_path, spacing_value);  break;
			case _USHORT: Linear_UpSample<unsigned short>(ref_path, input_path, output_path, spacing_value);  break;
			case _SHORT: Linear_UpSample<short>(ref_path, input_path, output_path, spacing_value);  break;
			case _UINT: Linear_UpSample<unsigned int>(ref_path, input_path, output_path, spacing_value);  break;
			case _INT: Linear_UpSample<int>(ref_path, input_path, output_path, spacing_value);  break;
			case _ULONG: Linear_UpSample<unsigned long>(ref_path, input_path, output_path, spacing_value);  break;
			case _LONG: Linear_UpSample<long>(ref_path, input_path, output_path, spacing_value);  break;
			case _FLOAT: Linear_UpSample<float>(ref_path, input_path, output_path, spacing_value);  break;
			case _DOUBLE: Linear_UpSample<double>(ref_path, input_path, output_path, spacing_value);  break;
			default: std::cout << "unrecognized pixel type"; break;
			}

		}
		else if (vm.count("upsample_self")) // upsample by itself
		{
			std::string input_path = Get<std::string>(vm, "in");
			std::string output_path = Get<std::string>(vm, "out");
			std::string spacing_file = Get<std::string>(vm, "spacing_file");

			std::cout << "input_path = " << input_path << std::endl;
			std::cout << "output_path = " << output_path << std::endl;
			std::cout << "spacing_file = " << spacing_file << std::endl;

			//UpSample_Self(input_path, output_path, spacing_file);

			LoadPixelType type = ReadPixelType(input_path.c_str());

			switch (type)
			{
			case _UCHAR: UpSample_Self<unsigned char>(input_path, output_path, spacing_file);  break;
			case _CHAR: UpSample_Self<char>(input_path, output_path, spacing_file);  break;
			case _USHORT: UpSample_Self<unsigned short>(input_path, output_path, spacing_file);  break;
			case _SHORT: UpSample_Self<short>(input_path, output_path, spacing_file);  break;
			case _UINT: UpSample_Self<unsigned int>(input_path, output_path, spacing_file);  break;
			case _INT: UpSample_Self<int>(input_path, output_path, spacing_file);  break;
			case _ULONG: UpSample_Self<unsigned long>(input_path, output_path, spacing_file);  break;
			case _LONG: UpSample_Self<long>(input_path, output_path, spacing_file);  break;
			case _FLOAT: UpSample_Self<float>(input_path, output_path, spacing_file);  break;
			case _DOUBLE: UpSample_Self<double>(input_path, output_path, spacing_file);  break;
			default: std::cout << "unrecognized pixel type"; break;
			}
		}
		else if (vm.count("label_upsample_self")) // upsample by label interpolation
		{
			std::string input_path = Get<std::string>(vm, "in");
			std::string output_path = Get<std::string>(vm, "out");
			std::string spacing_file = Get<std::string>(vm, "spacing_file");

			std::cout << "input_path = " << input_path << std::endl;
			std::cout << "output_path = " << output_path << std::endl;
			std::cout << "spacing_file = " << spacing_file << std::endl;

			Label_UpSample_Self(input_path, output_path, spacing_file);
		}
		else if (vm.count("resample")) // resample
		{
			std::string input_path = Get<std::string>(vm, "in");
			std::string output_path = Get<std::string>(vm, "out");

			std::cout << "input_path = " << input_path << std::endl;
			std::cout << "output_path = " << output_path << std::endl;

			Resample(input_path, output_path);
		}
		else if (vm.count("downsample")) //downsample
		{
			std::string input_path = Get<std::string>(vm, "in");
			double spacing_value = Get<double>(vm, "spacing");
			std::string output_path = Get<std::string>(vm, "out");

			std::cout << "input_path = " << input_path << std::endl;
			std::cout << "spacing_value = " << spacing_value << std::endl;
			std::cout << "output_path = " << output_path << std::endl;

			DownSample(input_path, spacing_value, output_path);
		}
		else if (vm.count("mask")) //mask
		{
			std::string input_path = Get<std::string>(vm, "in");
			std::string output_path = Get<std::string>(vm, "out");
			double threshold = Get<double>(vm, "threshold");
			
			std::cout << "input_path = " << input_path << std::endl;
			std::cout << "output_path = " << output_path << std::endl;

			GenegrateMaskImage(input_path, threshold, output_path);
		}
		else if (vm.count("world2voxel")) //world2voxel
		{
			std::string input_path = Get<std::string>(vm, "in");
			std::string landmark_path = Get<std::string>(vm, "landmark");

			std::cout << "input_path = " << input_path << std::endl;
			std::cout << "landmark = " << landmark_path << std::endl;

			World2Voxel(input_path, landmark_path);
		}
		else if (vm.count("voxel2world")) //voxel2world
		{
			std::string input_path = Get<std::string>(vm, "in");
			std::string landmark_path = Get<std::string>(vm, "landmark");

			std::cout << "input_path = " << input_path << std::endl;
			std::cout << "landmark = " << landmark_path << std::endl;

			Voxel2World(input_path, landmark_path);
		}
		else if (vm.count("histogram_match")) //histogram_match
		{
			std::string input_path = Get<std::string>(vm, "in");
			std::string ref_path = Get<std::string>(vm, "ref");
			std::string output_path = Get<std::string>(vm, "out");

			std::cout << "input_path = " << input_path << std::endl;
			std::cout << "ref_path = " << ref_path << std::endl;
			std::cout << "output_path = " << output_path << std::endl;

			HistogramMatching(ref_path, input_path, output_path);
		}
		else if (vm.count("generate_segment")) //generate_segment
		{
			std::string input_path = Get<std::string>(vm, "in");
			std::string output_path = Get<std::string>(vm, "out");
			double threshold = Get<double>(vm, "threshold");

			std::cout << "input_path = " << input_path << std::endl;
			std::cout << "output_path = " << output_path << std::endl;
			std::cout << "threshold = " << threshold << std::endl;
			
			GenegrateSegmentation(input_path, threshold, output_path);
		}
		else if (vm.count("largest_component")) //largest_component
		{
			std::string input_path = Get<std::string>(vm, "in");
			std::string output_path = Get<std::string>(vm, "out");

			std::cout << "input_path = " << input_path << std::endl;
			std::cout << "output_path = " << output_path << std::endl;

			ExtractLargestComponent(input_path, output_path);
		}
		else if (vm.count("fill_hole")) //fill_hole
		{
			std::string input_path = Get<std::string>(vm, "in");
			std::string output_path = Get<std::string>(vm, "out");
			int radius = Get<int>(vm, "radius");

			std::cout << "input_path = " << input_path << std::endl;
			std::cout << "output_path = " << output_path << std::endl;
			std::cout << "radius = " << radius << std::endl;

			HillHole(input_path, output_path, radius);

		}
		else if (vm.count("correct")) //correct
		{
			std::string input_path = Get<std::string>(vm, "in");
			std::string ref_path = Get<std::string>(vm, "ref");
			std::string output_path = Get<std::string>(vm, "out");
			double spacing_value = Get<double>(vm, "spacing");

			std::cout << "input_path = " << input_path << std::endl;
			std::cout << "ref_path = " << ref_path << std::endl;
			std::cout << "output_path = " << output_path << std::endl;
			std::cout << "spacing_value = " << spacing_value << std::endl;

			Correct(ref_path, input_path, output_path, spacing_value);
		}
		else if (vm.count("dilate")) //dilate
		{
			std::string input_path = Get<std::string>(vm, "in");
			std::string output_path = Get<std::string>(vm, "out");
			int radius = Get<int>(vm, "radius");

			std::cout << "input_path = " << input_path << std::endl;
			std::cout << "output_path = " << output_path << std::endl;
			std::cout << "radius = " << radius << std::endl;

			Dilate(input_path, output_path, radius);

		}
		else if (vm.count("erode")) //erode
		{
			std::string input_path = Get<std::string>(vm, "in");
			std::string output_path = Get<std::string>(vm, "out");
			int radius = Get<int>(vm, "radius");

			std::cout << "input_path = " << input_path << std::endl;
			std::cout << "output_path = " << output_path << std::endl;
			std::cout << "radius = " << radius << std::endl;

			Erode(input_path, output_path, radius);

		}
		else if (vm.count("convert")) //convert
		{
			std::string input_path = Get<std::string>(vm, "in");
			std::string output_path = Get<std::string>(vm, "out");

			std::cout << "input_path = " << input_path << std::endl;
			std::cout << "output_path = " << output_path << std::endl;

			Convert(input_path, output_path);

		}
		else if (vm.count("count")) //count
		{
			std::string input_path = Get<std::string>(vm, "in");

			std::cout << "input_path = " << input_path << std::endl;

			Count(input_path);

		}
		else if (vm.count("boundingbox")) //boundingbox
		{
			std::string input_path = Get<std::string>(vm, "in");
			std::string start_pt = Get<std::string>(vm, "start_pt");
			std::string end_pt = Get<std::string>(vm, "end_pt");

			std::cout << "input_path = " << input_path << std::endl;
			std::cout << "start_pt = " << start_pt << std::endl;
			std::cout << "end_pt = " << end_pt << std::endl;

			BoundingBox(input_path, start_pt, end_pt);
		}
		else if (vm.count("crop")) //crop
		{
			std::string input_path = Get<std::string>(vm, "in");
			std::string output_path = Get<std::string>(vm, "out");
			std::string start_pt = Get<std::string>(vm, "start_pt");
			std::string end_pt = Get<std::string>(vm, "end_pt");
			int pad = Get<int>(vm, "pad");

			std::cout << "input_path = " << input_path << std::endl;
			std::cout << "output_path = " << output_path << std::endl;
			std::cout << "start_pt = " << start_pt << std::endl;
			std::cout << "end_pt = " << end_pt << std::endl;
			std::cout << "pad = " << pad << std::endl;

			Crop(input_path, start_pt, end_pt, pad, output_path);

		}
		else if (vm.count("crop_femur")) //crop for femur
		{
			std::string input_path = Get<std::string>(vm, "in");
			std::string output_path = Get<std::string>(vm, "out");
			std::string landmark_path = Get<std::string>(vm, "landmark");

			std::cout << "input_path = " << input_path << std::endl;
			std::cout << "output_path = " << output_path << std::endl;
			std::cout << "landmark_path = " << landmark_path << std::endl;

			CropFemur(input_path, landmark_path, output_path);

		}
		else if (vm.count("minus")) //minus
		{
			std::string input_path = Get<std::string>(vm, "in");
			std::string output_path = Get<std::string>(vm, "out");
			int minus_value = Get<int>(vm, "minus_value");

			std::cout << "input_path = " << input_path << std::endl;
			std::cout << "minus_value = " << minus_value << std::endl;
			std::cout << "output_path = " << output_path << std::endl;

			Minus(input_path, minus_value, output_path);
		}
		else if (vm.count("putback")) //putback
		{
			std::string ref_path = Get<std::string>(vm, "ref");
			std::string input_path = Get<std::string>(vm, "in");
			std::string output_path = Get<std::string>(vm, "out");
			std::string start_pt = Get<std::string>(vm, "start_pt");
			std::string end_pt = Get<std::string>(vm, "end_pt");
			
			std::cout << "ref_path = " << input_path << std::endl;
			std::cout << "input_path = " << input_path << std::endl;
			std::cout << "output_path = " << output_path << std::endl;
			std::cout << "start_pt = " << start_pt << std::endl;
			std::cout << "end_pt = " << end_pt << std::endl;

			PutBack(ref_path, input_path, start_pt, end_pt, output_path);
		}
		else if (vm.count("label_transform")) //label_transform
		{
			std::string input_path = Get<std::string>(vm, "in");
			std::string output_path = Get<std::string>(vm, "out");

			std::cout << "input_path = " << input_path << std::endl;
			std::cout << "output_path = " << output_path << std::endl;

			LabelTransform(input_path, output_path);
		}
		else if (vm.count("adaptive_histogram_equalization")) //adaptive_histogram_equalization
		{
			std::string input_path = Get<std::string>(vm, "in");
			std::string output_path = Get<std::string>(vm, "out");
			int radius = Get<int>(vm, "radius");
			double alpha = Get<double>(vm, "alpha");
			double beta = Get<double>(vm, "beta");

			std::cout << "input_path = " << input_path << std::endl;
			std::cout << "output_path = " << output_path << std::endl;
			std::cout << "radius = " << radius << std::endl;
			std::cout << "alpha = " << alpha << std::endl;
			std::cout << "beta = " << beta << std::endl;


			AdaptiveHistogramEqualization(input_path, output_path, alpha, beta, radius);
		}
		else if (vm.count("rescale_intensity")) //rescale_intensity
		{
			std::string input_path = Get<std::string>(vm, "in");
			std::string output_path = Get<std::string>(vm, "out");
			int min = Get<int>(vm, "min");
			int max = Get<int>(vm, "max");

			std::cout << "input_path = " << input_path << std::endl;
			std::cout << "output_path = " << output_path << std::endl;
			std::cout << "min = " << min << std::endl;
			std::cout << "max = " << max << std::endl;


			RescaleIntensity(input_path, min, max, output_path);
		}
		else if (vm.count("orient")) //orient
		{
			std::string input_path = Get<std::string>(vm, "in");
			std::string output_path = Get<std::string>(vm, "out");

			std::cout << "input_path = " << input_path << std::endl;
			std::cout << "output_path = " << output_path << std::endl;

			Orient(input_path, output_path);
		}
		else if (vm.count("judge")) //judge
		{
			std::string input_path = Get<std::string>(vm, "in");		
			std::string ref_path = Get<std::string>(vm, "ref");

			std::cout << "input_path = " << input_path << std::endl;
			std::cout << "ref_path = " << ref_path << std::endl;

			Judge(ref_path, input_path);
		}
		else if (vm.count("cutoff"))//cutoff
		{
			std::string input_path = Get<std::string>(vm, "in");
			std::string output_path = Get<std::string>(vm, "out");
			double ratio = Get<double>(vm, "ratio");

			std::cout << "input_path = " << input_path << std::endl;
			std::cout << "output_path = " << output_path << std::endl;
			std::cout << "ratio = " << ratio << std::endl;

			CutOff(input_path, output_path, ratio);
		}
		else if (vm.count("extract_landmark_femur"))//extract_landmark_femur
		{
			std::string input_path = Get<std::string>(vm, "in");
			std::string output_path = Get<std::string>(vm, "out");
			std::string start_pt = Get<std::string>(vm, "start_pt");
			std::string end_pt = Get<std::string>(vm, "end_pt");

			std::cout << "input_path = " << input_path << std::endl;
			std::cout << "output_path = " << output_path << std::endl;
			std::cout << "start_pt = " << start_pt << std::endl;
			std::cout << "end_pt = " << end_pt << std::endl;

			ExtractLandmarkOfFemur(input_path, start_pt, end_pt, output_path);
		}
		else
		{
			std::cerr << "option not supported ( must pick one option in functional options list )" << std::endl;
			return -1;
		}
		
	}
	catch (boost::program_options::error& exp) {
		std::cerr << exp.what() << std::endl;
		return -1;
	}

return 0;
}
