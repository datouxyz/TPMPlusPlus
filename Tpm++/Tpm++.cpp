// Tpm++.cpp: 定义应用程序的入口点。
//

#include "Tpm++.h"
#include "TissueProblityMap.h"
int GetNumPhyCpu();
#include "TPMInternalFuncs.h"
using namespace std;
//#include "struct2x/json/decoder.h"
//#include "struct2x/json/encoder.h"
#include <fstream>
// enum TPMOutMask
// {
// AffineMat=0x00000001,FieldInverse=0x00000002,FieldForward=0x00000004,BiasCorrected=0x00000008,
// BiasField=0x00000010,TissImported=0x00000020,TissNative =0x00000040,TissModulate=0x00000080,
// TissUnModulate =0x00000100
// };
std::vector<std::string> SplitString(const std::string& s, char delimiter)
{
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter))
    {
        tokens.push_back(token);
    }
    return tokens;
}
sTissues DefaultTissues()
{
	sTissues tss=
	{sTissue(1, true, false, false, false), sTissue(1, true, false, false, false),
		sTissue(2, true, false, false, false), sTissue(3, true, false, false, false),
		sTissue(4, true, false, false, false), sTissue(2, false, false, false, false) };
	return tss;
}

#include <sitkExceptionObject.h>
#include <ITK-prefix/include/ITK-5.2/itkMacro.h>

sitk::Transform GetTransformFromAffine(MAT_TYPE Affine)
{	
	auto transform =sitk::AffineTransform(3);	
		
	std::vector<double> parameters(12);
	//Copy the Eigen matrix to the ITK matrix and offset
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			parameters[i*3 + j] =Affine (i, j);
		}
		parameters[9 + i] = Affine(i, 3);
	}
	
	transform.SetParameters(parameters);

	return transform;
}

// "E:\x64\XCrossDebug\SimpleITK_SimpleITKFilters-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITKBasicFilters0-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITKBasicFilters1-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITKCommon-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITKIO-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITKRegistration-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKAnisotropicSmoothing-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKAntiAlias-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKBiasCorrection-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKBinaryMathematicalMorphology-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKClassifiers-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKColormap-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKCommon-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKConnectedComponents-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKConvolution-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKCurvatureFlow-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKDeconvolution-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKDenoising-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKDisplacementField-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKDistanceMap-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKFastMarching-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKFFT-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKImageCompare-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKImageCompose-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKImageFeature-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKImageFilterBase-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKImageFunction-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKImageFusion-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKImageGradient-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKImageGrid-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKImageIntensity-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKImageLabel-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKImageNoise-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKImageSources-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKImageStatistics-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKLabelMap-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKLabelVoting-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKLevelSets-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKMathematicalMorphology-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKPDEDeformableRegistration-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKRegionGrowing-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKRegistrationCommon-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKReview-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKSmoothing-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKSuperPixel-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKThresholding-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKTransform-2.1.dll"
// "E:\x64\XCrossDebug\SimpleITK_ITKWatersheds-2.1.dll"
// {
// SimpleITK_SimpleITKFilters,
// SimpleITKBasicFilters0-2.1,
// SimpleITKBasicFilters1-2.1,


// }
inline void LoadVecFromArray(JSONNode & arraynode, FloatVec&outvec)
{
	outvec.clear();
	for (auto & el : arraynode)
	{
		outvec.push_back(el.as_float());
	}
}
inline void LoadVecFromArray(JSONNode & arraynode, DoubleVec&outvec)
{
	outvec.clear();
	for (auto & el : arraynode)
	{
		outvec.push_back(el.as_float());
	}
}
inline void SaveVecToArray(JSONNode & arraynode, FloatVec&invec)
{
	for (auto val : invec)
	{
		arraynode.push_back(JSONNode("", val));
	}
}
inline void SaveVecToArray(JSONNode & arraynode, DoubleVec&invec)
{
	for (auto val : invec)
	{
		arraynode.push_back(JSONNode("", val));
	}
}
typedef std::vector<string> StringVec;
inline void SaveVecToArray(JSONNode & arraynode, StringVec &invec)
{
	for (auto val : invec)
	{
		arraynode.push_back(JSONNode("", val));
	}
}


void LoadTPMInputFromJson(string & strjson, sTpmInputs &inputs)
{
	JSONNode root = libjson::parse(strjson.c_str());
	auto & tisses = root["m_tissues"].as_array();
	inputs.m_tissues.clear();
	for (auto &tiss : tisses)
	{
		sTissue newtis;
		newtis.m_ngaus = tiss["m_ngaus"].as_int();
		newtis.m_NativeTissue_NativeSpace = tiss["m_NativeTissue_NativeSpace"].as_bool();
		newtis.m_NativeTissue_DartelImported = tiss["m_NativeTissue_DartelImported"].as_bool();
		newtis.m_WarpTissue_Modulated = tiss["m_WarpTissue_Modulated"].as_bool();
		newtis.m_WarpTissue_UnModulated = tiss["m_WarpTissue_UnModulated"].as_bool();
		inputs.m_tissues.push_back(newtis);
	}
	auto & channelnodes = root["m_Channels"].as_array();
	inputs.m_Channels.clear();
	for (auto &chann : channelnodes)
	{
		sChannel newchan;
		newchan.m_Bisareg = chann["m_Bisareg"].as_float();
		newchan.m_Biasfwhm = chann["m_Biasfwhm"].as_float();
		newchan.m_biascorrect = chann["m_biascorrect"].as_bool();
		newchan.m_biasfield = chann["m_biasfield"].as_float();
		inputs.m_Channels.push_back(newchan);
	}

	inputs.m_nsample = root["m_nsample"].as_int();
	inputs.m_fwhm = root["m_fwhm"].as_float();

	LoadVecFromArray(root["m_reg"].as_array(), inputs.m_reg);
	LoadVecFromArray(root["m_lkp"].as_array(), inputs.m_lkp);
	inputs.m_ptype = (PriorsType)root["m_ptype"].as_int();
	inputs.m_dFieldInverse = root["m_dFieldInverse"].as_bool();
	inputs.m_dFieldForward = root["m_dFieldForward"].as_bool();
	inputs.m_Mrf = root["m_Mrf"].as_bool();
	inputs.m_Cleanup = root["m_Cleanup"].as_int();
	inputs.m_dFieldNative = root["m_dFieldNative"].as_int();

}


void SaveTPMInputsToJson(string & strjson, sTpmInputs &inputs)
{
	JSONNode root;

	{
		JSONNode tissnodes(JSON_ARRAY);	tissnodes.set_name("m_tissues");
		for (auto & tiss : inputs.m_tissues)
		{
			JSONNode ntiss;

			ntiss.push_back(JSONNode("m_ngaus", tiss.m_ngaus));
			ntiss.push_back(JSONNode("m_NativeTissue_NativeSpace", tiss.m_NativeTissue_NativeSpace));
			ntiss.push_back(JSONNode("m_NativeTissue_DartelImported", tiss.m_NativeTissue_DartelImported));
			ntiss.push_back(JSONNode("m_WarpTissue_Modulated", tiss.m_WarpTissue_Modulated));
			ntiss.push_back(JSONNode("m_WarpTissue_UnModulated", tiss.m_WarpTissue_UnModulated));
			tissnodes.push_back(ntiss);
		}
		root.push_back(tissnodes);
	}

	{
		JSONNode channelnodes(JSON_ARRAY);	channelnodes.set_name("m_Channels");
		for (auto & chan : inputs.m_Channels)
		{
			JSONNode nch;
			nch.push_back(JSONNode("m_Bisareg", chan.m_Bisareg));
			nch.push_back(JSONNode("m_Biasfwhm", chan.m_Biasfwhm));
			nch.push_back(JSONNode("m_biascorrect", chan.m_biascorrect));
			nch.push_back(JSONNode("m_biasfield", chan.m_biasfield));
			channelnodes.push_back(nch);
		}
		root.push_back(channelnodes);
	}

	root.push_back(JSONNode("m_nsample", inputs.m_nsample));
	root.push_back(JSONNode("m_fwhm", inputs.m_fwhm));

	{
		JSONNode regarray(JSON_ARRAY); regarray.set_name("m_reg");
		SaveVecToArray(regarray, inputs.m_reg);
		root.push_back(regarray);
	}
	{
		JSONNode lkparray(JSON_ARRAY); lkparray.set_name("m_lkp");
		SaveVecToArray(lkparray, inputs.m_lkp);
		root.push_back(lkparray);
	}

	root.push_back(JSONNode("m_ptype", inputs.m_ptype));
	root.push_back(JSONNode("m_dFieldInverse", inputs.m_dFieldInverse));
	root.push_back(JSONNode("m_dFieldForward", inputs.m_dFieldForward));
	root.push_back(JSONNode("m_Mrf", inputs.m_Mrf));
	root.push_back(JSONNode("m_Cleanup", inputs.m_Cleanup));
	root.push_back(JSONNode("m_dFieldNative", inputs.m_dFieldNative));

	strjson = root.write();
}




int main(int argc, char *argv[])
{
	if (argc < 5)
	{
		cout << "Usage: " << argv[0] << " -in inputimage -tpm tpmimage -out outFile -json TpmParam " << endl;		
		 cout<<"Out file will has a postfix 'outimage_FieldInverse.nii, outimage_AffineMat.txt,etc...'"<<endl;
		return 1;
	}
	sChannels Channels(1);
	sTpmInputs  TpmInputs;
	ItkImageVec InputImages;
	sitk::Image TpmImg;
	string outimgpath;
	bool bAffineOnly=false;
	// check params for -in -tpm -out -outaffine in a loop
	for(int iarg = 1;iarg<argc;iarg++)
	{
		if (strcmp(argv[iarg], "-in") == 0)
		{
			if (iarg + 1 >= argc)
			{
				cout << "missing input image" << endl;
				return 1;
			}
			// load images until next '-'param
			while (iarg + 1 < argc && argv[iarg + 1][0] != '-')
			{
				sitk::Image img = sitk::ReadImage(argv[iarg+1]);
				InputImages.push_back(img);
				iarg++;
			}
		}
		else if (strcmp(argv[iarg], "-tpm") == 0)
		{
			if (iarg + 1 >= argc)
			{
				cout << "missing tpm image" << endl;
				return 1;
			}
			TpmImg = sitk::ReadImage(argv[iarg+1],sitk::sitkFloat64);
			iarg++;
		}
		else if (strcmp(argv[iarg], "-out") == 0)
		{
			if (iarg + 1 >= argc)
			{
				cout << "missing out outPath" << endl;
				return 1;
			}
			outimgpath = argv[iarg+1];
			iarg++;
		}
		else if (strcmp(argv[iarg], "-json") == 0)
		{
			if (iarg + 1 >= argc)
			{
				cout << "missing json params "<< endl;
				return 1;
			}		
			string jsonstr = LoadStringFromFile(argv[iarg+1]);			

			LoadTPMInputFromJson(jsonstr,TpmInputs);
			//bool bDecode = struct2x::JSONDecoder(jsonstr.c_str()) >> TpmInputs;			
			iarg++;
		}
		else if (strcmp(argv[iarg], "-AffineOnly") == 0)
		{
			if (iarg + 1 >= argc)
			{
				cout << "missing AffineOnly value"<< endl;
				bAffineOnly = false;
			}
			else
			bAffineOnly = (argv[iarg+1] == string("1")) ? true:false;
		}
	}
	
	if(InputImages.size()==0)
	{	
		cout << "missing input image" << endl;
		return 1;
	}
	if(TpmImg.GetSize().size()==0)
	{
		cout << "missing tpm image" << endl;
		return 1;
	}
	if(outimgpath.size()==0)
	{
		cout << "missing out outimgpath" << endl;
		return 1;
	}	

	// {
	// 	MAT_TYPE Affine(4,4);
	// 	Affine<< 1.05645,    -0.0205341, -0.0266813, 0.563406,
	// 	0.0177077,  0.998854,   0.305593,  2.91475,
	// 	0.0116576,  -0.350093,  1.08758,   -25.5842,
	// 	0,          0,          0,         1;

	// 	MAT_TYPE flipMatrix = MAT_TYPE::Identity(4,4);
	// 	flipMatrix(0, 0) = -1;
	// 	flipMatrix(1, 1) = -1;
	// 	Affine = flipMatrix * Affine * flipMatrix;
	// 	cout<<Affine<<endl;
	//     auto & movimg =  InputImages[0];
	// 	sitk::WriteImage(sitk::Resample(movimg,GetTransformFromAffine(Affine.inverse())),"e:\\movimg.nii" );
	// }

	Eigen::initParallel();
	int ncpu = GetNumPhyCpu();
	Eigen::setNbThreads(ncpu);
	omp_set_num_threads(ncpu);
	while (InputImages.size() > TpmInputs.m_Channels.size())
	{
		TpmInputs.m_Channels.push_back(sChannel());
	}
	if (TpmImg.GetPixelID() != sitk::sitkFloat64)
	{
		ErrOutput(u8"Warn",u8"TPM Will be cast to Float 64(8byte)");
		TpmImg = sitk::Cast(TpmImg,sitk::sitkFloat64);		
	}

	sTPMPriors tpmPriors;
	tpmPriors.LoadPriors(TpmImg);
	
	for (int ichan = 0; ichan < TpmInputs.m_Channels.size(); ichan++)
	{
		auto &chan = TpmInputs.m_Channels[ichan];
		auto & inimg = InputImages[ichan];
		if (inimg.GetPixelID() != sitk::sitkFloat32)
		{
			ErrOutput( u8"Warn",u8"Input Image will be cast to Float 32(4byte)");
			inimg = sitk::Cast(inimg,sitk::sitkFloat32);
		}
		chan.m_Img = inimg;
	}
	
	MAT_TYPE Affine = TpmInputs.m_InitAffine;
	cout<<Affine<<endl;
	ErrOutput(u8"Affine:", u8"AffinePass Start");
	
	Affine = AffinePass(TpmInputs.m_Channels[0].m_Img, TpmInputs.m_nsample, (TpmInputs.m_fwhm + 1) * 16, tpmPriors, Affine, TpmInputs.m_ptype);
	cout<<Affine<<endl;
	Affine = AffinePass(TpmInputs.m_Channels[0].m_Img, TpmInputs.m_nsample, TpmInputs.m_fwhm, tpmPriors, Affine, TpmInputs.m_ptype);
	cout<<Affine<<endl;
	
	ErrOutput(u8"Affine:", u8"AffinePass Fin");
 
	if(outimgpath.size()>0)
	{
		MAT_TYPE flipMatrix = MAT_TYPE::Identity(4,4);
		flipMatrix(0, 0) = -1;
		flipMatrix(1, 1) = -1;
		auto ITKAffine = flipMatrix * Affine * flipMatrix;

		ofstream affinefile(outimgpath + "_AffineMat.txt");
		cout << ITKAffine<<endl;;
		affinefile << ITKAffine;			
	}

	if(bAffineOnly)
	{
		ErrOutput(u8"Affine:","Affine Only Done quit;");
		return 0;
	}

	sWrapResults wpres;
	WrapPass(TpmInputs, tpmPriors, Affine, wpres);

	sWriteImgs sWrites;
	TpmWrite8(TpmInputs, tpmPriors, wpres, sWrites);

	if (!sWrites.FieldInverse.IsEmpty() )
	{
		sitk::WriteImage(sWrites.FieldInverse.m_Img, outimgpath + "_FieldInverse.nii");
		cout<<outimgpath + "_FieldInverse.nii"<<" Writed"<<endl;
	}
	//if (std::find(modes.begin(),modes.end(),"FieldForwared")!=modes.end())
	if(!sWrites.FieldForwared.IsEmpty())
	{
		sitk::WriteImage(sWrites.FieldForwared.m_Img, outimgpath + "_FieldForwared.nii");
		cout<<outimgpath + "_FieldForwared.nii"<<" Writed"<<endl;
	}
	
	//if (std::find(modes.begin(),modes.end(),"BiasCorrected")!=modes.end())
	//if(sWrites.)	
	for (int ichan=0;ichan<sWrites.chans.size();ichan++)
	{
		auto & chan = sWrites.chans[ichan];
		if (!chan.m_Nc.IsEmpty())
		{
			auto outfile = FormatedString("%s_BiasCorrected_%d.nii",outimgpath.c_str(),ichan);
			sitk::WriteImage(chan.m_Nc.m_Img, outfile);
			cout<<outfile<<" Writed"<<endl;
		}
	}
	{
		for (int ichan=0;ichan<sWrites.chans.size();ichan++)
		{
			auto & chan = sWrites.chans[ichan];
			if (!chan.m_Nf.IsEmpty())
			{
				auto outfile = FormatedString("%s_BiasField_%d.nii",outimgpath.c_str(),ichan);
				sitk::WriteImage(chan.m_Nf.m_Img, outfile);				
				cout<<outfile<<" Writed"<<endl;
			}
		}		
	}	
	{
		for (int ichan=0;ichan<sWrites.tissImported.size();ichan++)
		{
			auto & ti = sWrites.tissImported[ichan];
			if (!ti.IsEmpty())
			{
				auto outfile = FormatedString("%s_TissImported_%d.nii",outimgpath.c_str(),ichan);
				sitk::WriteImage(ti.m_Img, outfile);
				cout<<outfile<<" Writed"<<endl;
			}
		}		
	}	
	{
		for (int ichan=0;ichan<sWrites.tissNative.size();ichan++)
		{
			auto & tn = sWrites.tissNative[ichan];
			if (!tn.IsEmpty())
			{
				//imgs.push_back(chan.m_Img);
				sitk::WriteImage(tn.m_Img, FormatedString("%s_TissNative_%d.nii",outimgpath.c_str(),ichan));
			}
		}		
	}

	//if (std::find(modes.begin(),modes.end(),"TissModulate")!=modes.end())
	{
		for (int ichan=0;ichan<sWrites.tissModulate.size();ichan++)
		{
			auto & tm = sWrites.tissModulate[ichan];
			if (!tm.IsEmpty())
			{
				//imgs.push_back(chan.m_Img);
				sitk::WriteImage(tm.m_Img, FormatedString("%s_TissModulate_%d.nii",outimgpath.c_str(),ichan));
			}
		}		
	}
	//if (std::find(modes.begin(),modes.end(),"TissUnModulate")!=modes.end())
	{
		for (int ichan=0;ichan<sWrites.tissUnModulate.size();ichan++)
		{
			auto & tu = sWrites.tissUnModulate[ichan];
			if (!tu.IsEmpty())
			{
				//imgs.push_back(chan.m_Img);
				sitk::WriteImage(tu.m_Img, FormatedString("%s_TissUnModulate_%d.nii",outimgpath.c_str(),ichan));
			}
		}		
	}	
	return 0;
}
// cout<<"Affine.inverse()\n"<<Affine.inverse()<<endl;
// 		cout<<"MG1.inverse()\n"<<MG1.inverse()<<endl;
// 		cout<<"MG1()\n"<<MG1<<endl;
// 		cout<<"Affine.inverse()*MG1.inverse()\n"<<Affine.inverse()*MG1.inverse()<<endl;
// 		cout<<" MG1* Affine.inverse()\n"<< MG1* Affine.inverse()<<endl;
// 		cout<<"Affine"<<Affine<<endl; 
// 		cout<<"MG1"<<MG1<<endl; 
	
// 		cout<<"TPMMat"<<TPMMat<<endl;
// //auto MG = GetMatFromImageDirOriSpace(movimg.GetDirection(), movimg.GetOrigin(), movimg.GetSpacing(), true, true);
	
		// //Affine= MG.inverse()*Affine.inverse();//(tpmPriors.m_TPMMat.inverse()*Affine*MG).inverse();
		// //convert Affine to ITK format
		// //Affine = Affine.inverse();
		// //Affine = MG*Affine*tpmPriors.m_TPMMat;
		// //Affine = Affine.inverse();
			
		// auto & movimg =  TpmInputs.m_Channels[0].m_Img;
		
		// auto size = TpmImg.GetSize();

		// std::vector<sitk::Image> Tpmimages;
		
		// sitk::Image TPM0_ = sitk::Extract(TpmImg, {size[0],size[1],size[2],0}, {0,0,0,0});
		// sitk::Image TPM0 = sitk::Cast(TPM0_,movimg.GetPixelID());

		// //movimg = sitk::Resample(movimg, TPM0);
		// //sitk::WriteImage(image3D,"e:\\outtest.nii");		
		// sitk::Transform transform=sitk::AffineTransform(3);
		//transform =  sitk::CenteredTransformInitializer(movimg,TPM0,sitk::AffineTransform(3));		
		// catch(const std::exception& e)
		// {
		// 	std::cerr << e.what() << '\n';
		// }
		// catch (itk::ExceptionObject & exp)
		// {
		// 	cout<< exp.GetDescription();			
		// }
		// catch (sitk::GenericException & exp)
		// {		
		// 	cout<<exp.GetDescription();
		// }
		// catch(...)
		// {			
		// 	return -1;
		// }
	
		//sitk::Image newimg = sitk::Resample(movimg, sitk::Transform(transform), sitk::sitkLinear, 0.0, movimg.GetPixelID());
		//sitk::Transform identityTransform = sitk::Transform();  // Identity transform
		//sitk::Image referenceImage({400,400,400},sitk::sitkFloat32);  // Copy the original image
		// referenceImage.SetDirection({1, 0, 0, 0, 1, 0, 0, 0, 1});  // Set direction to (1,1,1)
		// referenceImage.SetSpacing({1, 1, 1});  // Set spacing to (1,1,1)
		// referenceImage.SetOrigin({-150, -150, -150});  // Set origin to (0,0,0)
			//cout<<"MG"<<MG<<endl; 
		//cout<<"MG1"<<MG1<<endl; 
		
		//auto transform =  sitk::CenteredTransformInitializer(movimg,TPM0,sitk::AffineTransform(3));	
			// vector<double> center = movimg.TransformContinuousIndexToPhysicalPoint(
		// 	{movimg.GetSize()[0] * 0.5,movimg.GetSize()[1] * 0.5,movimg.GetSize()[2] * 0.5}
		// );
		// transform.SetCenter(center);
	

		//sitk::Image newimg = sitk::Resample(movimg, referenceImage,transform);
// //auto ApplyTr = (TPMMat * Affine.inverse()*MG.inverse()).inverse();
		// auto ApplyTr = (MG.inverse() * Affine*MG).inverse();
		//auto ApplyTr = TPMMat.inverse()*Affine*MG;
		//auto ApplyTr = (TPMMat.inverse()*Affine*MG).inverse();//*Affine*TPMMat.inverse();//TPMMat1.inverse()*Affine*MG1;
		//MG.inverse()*Affine.inverse()*TPMMat; //
		//auto ApplyTr = (TPMMat1.inverse()*Affine*MG1).inverse();
		//auto ApplyTr = (TPMMat1.inverse()*Affine*MG1);
		//auto ApplyTr = (MG1.inverse() *Affine.inverse()*TPMMat1);
		//auto ApplyTr = (MG1.inverse() *Affine*TPMMat1);
		//auto ApplyTr = (Affine*TPMMat1);
		//auto ApplyTr = (Affine*TPMMat1.inverse());
		//auto ApplyTr = (MG1.inverse() * Affine*MG1);
		//auto ApplyTr = (MG1 * Affine*MG1.inverse());
		//auto ApplyTr = (MG1.inverse() * Affine*MG1);
		//auto ApplyTr = (MG1.inverse() * Affine.inverse()*MG1);		
		//T = priors.m_TPMMat.inverse()*M*regdata.m_MG;
		// Affine(0,3)=1;
		// Affine(1,3)=1;
		// Affine(2,3)=1;
		//auto ApplyTr = (MG1.inverse()* Affine*MG1).inverse(); /// angle matched offset prob
		//auto ApplyTr = TPMMat* Affine.inverse()*TPMMat.inverse(); /// angle matched offset prob
		
		//auto ApplyTr = (TPMMat.inverse()*Affine*MG1).inverse();// nop
		//auto ApplyTr = MG1.inverse()*Affine.inverse()*TPMMat;// nop
		//auto ApplyTr = MG* Affine.inverse()*MG.inverse();
		//auto ApplyTr = MG.inverse()* Affine.inverse()*MG;
		//auto ApplyTr = MG.inverse()* Affine*MG;
		//auto ApplyTr =MG1*Affine.inverse()*MG1.inverse(); 
		//auto ApplyTr =MG1.inverse()*Affine.inverse()*MG1; //offset prob
		//auto ApplyTr = TPMMat.inverse()*Affine.inverse()*TPMMat; //offset prob
		//auto ApplyTr = MG1*Affine.inverse()*TPMMat.inverse(); 
		//auto ApplyTr = MG1.inverse()*Affine.inverse()*TPMMat; 
		//auto ApplyTr = MG1.inverse()* Affine.inverse()*MG1; ///// angle matched wrong scale
		//auto ApplyTr = TPMMat1 * Affine.inverse()*MG1.inverse(); 
		//auto ApplyTr = TPMMat1* Affine.inverse()*TPMMat1.inverse(); 
		
		//auto ApplyTr = TPMMat1* Affine.inverse()*MG1.inverse(); 
		//auto ApplyTr = MG1* Affine.inverse()*TPMMat1.inverse();
		//auto ApplyTr = MG* Affine.inverse()*MG.inverse();
		//auto ApplyTr = MG.inverse()* Affine.inverse()*MG;
		//auto ApplyTr = TPMMat.inverse()* Affine.inverse()*MG;