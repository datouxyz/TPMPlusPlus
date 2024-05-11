# TPM++
TPM++ Read Me

TPM++ is a C++ rewrite of the Tissue Probability Map in SPM12. It uses Eigen for matrix operations to replace the Matlab framework and uses SimpleITK as the image library. The image format and coordinate space are consistent with ITK. TPM++ is based on the GPL v2 license. TPM++ produces the same results as SPM12, but its computational performance is N times faster than SPM12, where N is the number of physical CPU cores. 
By default, TPM++ does not use hyper-threading, as hyper-threading can no longer provide high performance.
TPM++ provides a command-line program to calculate Wrap and perform image segmentation: The command is as follows: 

**TPM++.exe -in image.nii -tpm D:/TPM.nii -out D:/_out.nii -json D:/tpminput.json -AffineOnly 0** <br>
Test Case:<br>
Input： <br>
![Alt Text](https://github.com/datouxyz/TPMPlusPlus/blob/master/ReadMe.files/ReadMe3193.png) <br>
Outputs:<br>

![Alt Text](https://github.com/datouxyz/TPMPlusPlus/blob/master/ReadMe.files/ReadMe3251.png) ![Alt Text](https://github.com/datouxyz/TPMPlusPlus/blob/master/ReadMe.files/ReadMe3253.png)
![Alt Text](https://github.com/datouxyz/TPMPlusPlus/blob/master/ReadMe.files/ReadMe3274.png) ![Alt Text](https://github.com/datouxyz/TPMPlusPlus/blob/master/ReadMe.files/ReadMe3276.png)<br>
Parameter description:<br>
-in: Input image path<br>
-tpm: TPM template path<br>
-out: Output file path, the actual output file will add a suffix to this path:<br>
                      _FieldInverse: Inverse deformation field<br>
                      _FieldForward: Forward deformation field<br>
                      _BiasCorrected: After bias correction<br>
                      _BiasField: Bias field<br>
                      _TissNative: Segmentation result (in original space)<br>
                      _TissImported: Segmentation result (in Dartel space)<br>
                      _TissModulate: Segmentation result (in imported space, after modulation)<br>
                      _TissUnModulate: Segmentation result (in imported space, before modulation)<br>
-AffineOnly: //When set to 1, only the affine transformation is calculated and the affine matrix is output.<br>
-json: The json file that records all parameters of TPM operation. <br>
The TpmInput.json file is a basic parameter case with a sample number of 3 and divided into 6 tissue types.<br>
Jsons sample：See TpmInput.json<br>
Jsons sample：<br>
{<br>
// "Settings for whether to generate related TPM for each tissue"<br>
"m_tissues": [<br>
{<br>
"m_ngaus": 1<br>
"m_NativeTissue_NativeSpace": true,<br>
"m_NativeTissue_DartelImported": true,<br>
"m_WarpTissue_Modulated": true,<br>
"m_WarpTissue_UnModulated": true<br>
},<br>
{<br>
"m_ngaus": 1,<br>
"m_NativeTissue_NativeSpace": true,<br>
"m_NativeTissue_DartelImported": true,<br>
"m_WarpTissue_Modulated": true,<br>
"m_WarpTissue_UnModulated": true<br>
},<br>
{<br>
"m_ngaus": 2,<br>
"m_NativeTissue_NativeSpace": true,<br>
"m_NativeTissue_DartelImported": true,<br>
"m_WarpTissue_Modulated": true,<br>
"m_WarpTissue_UnModulated": true<br>
},<br>
{<br>
"m_ngaus": 3,<br>
"m_NativeTissue_NativeSpace": true,<br>
"m_NativeTissue_DartelImported": true,<br>
"m_WarpTissue_Modulated": true,<br>
"m_WarpTissue_UnModulated": true<br>
},<br>
{<br>
"m_ngaus": 4,<br>
"m_NativeTissue_NativeSpace": true,<br>
"m_NativeTissue_DartelImported": true,<br>
"m_WarpTissue_Modulated": true,<br>
"m_WarpTissue_UnModulated": true<br>
},<br>
{<br>
"m_ngaus": 2,<br>
"m_NativeTissue_NativeSpace": true,<br>
"m_NativeTissue_DartelImported": true,<br>
"m_WarpTissue_Modulated": true,<br>
"m_WarpTissue_UnModulated": true<br>
}<br>
],<br>
"m_Channels": [// Settings for each input file (channel) <br>
{<br>
"m_Bisareg": 0.001000,<br>
"m_Biasfwhm": 60.000000,<br>
"m_biascorrect": true,<br>
"m_biasfield": true<br>
}<br>
],<br>
"m_nsample": 3,//sample count<br>
"m_fwhm": 0.000000,//smooth<br>
"m_reg": [//Regularization <br>
0.000000,<br>
0.001000,<br>
0.500000,<br>
0.050000,<br>
0.200000<br>
],<br>
"m_lkp": [//Gaussian function lookup table, sTpmInputs(sTissues &tss, sChannels &Channels) constructor can calculate this table  <br>
   1.000000,<br>
2.000000,<br>
3.000000,<br>
3.000000,<br>
4.000000,<br>
4.000000,<br>
4.000000,<br>
5.000000,<br>
5.000000,<br>
5.000000,<br>
5.000000,<br>
6.000000,<br>
6.000000<br>
],<br>
"m_ptype": 0,<br>
/*Prior type：<br>
enum PriorsType<br>
{<br>
PriorsType_mni = 0,<br>
PriorsType_imni,<br>
PriorsType_rigid,<br>
PriorsType_subject,<br>
PriorsType_eastern,<br>
PriorsType_none = 0xffffffff<br>
};*/<br>
"m_dFieldInverse": true,<br>
"m_dFieldForward": true,<br>
"m_Mrf": false,<br>
"m_Cleanup": 1<br>
}<br>
Folder description:<br>
Eigen-3.4.0: Eigen dependency package <br>
ITK-prefix: ITK header files <br>
Sitk: SITK header files <br>
Libjson: jsonlibjson header files <br>
SPMSRC: SPM native C code <br>
Struct2x: Struct2x serialization library



