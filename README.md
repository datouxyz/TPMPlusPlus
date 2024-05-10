# TPMPP
TPM++ Read Me

TPM++ is a C++ rewrite of the Tissue Probability Map in SPM12. It uses Eigen for matrix operations to replace the Matlab framework and uses SimpleITK as the image library. The image format and coordinate space are consistent with ITK. TPM++ is based on the GPL v2 license. TPM++ produces the same results as SPM12, but its computational performance is N times faster than SPM12, where N is the number of physical CPU cores. 
By default, TPM++ does not use hyper-threading, as hyper-threading can no longer provide high performance.
TPM++ provides a command-line program to calculate Wrap and perform image segmentation: The command is as follows: 

TPM++.exe -in image.nii -tpm D:/TPM.nii -out D:/_out.nii -json D:/tpminput.json -AffineOnly 0 

Parameter description:
-in: Input image path 
-tpm: TPM template path 
-out: Output file path, the actual output file will add a suffix to this path: 
           _FieldInverse: Inverse deformation field.
           
           _FieldForward: Forward deformation field.
           
           _BiasCorrected: After bias correction.
           
           _BiasField: Bias field.
           
           _TissNative: Segmentation result (in original space).
           
           _TissImported: Segmentation result (in Dartel space).
           
           _TissModulate: Segmentation result (in imported space, after modulation)
           
           _TissUnModulate: Segmentation result (in imported space, before modulation)
           
-AffineOnly: //When set to 1, only the affine transformation is calculated and the affine matrix is output.

-json: The json file that records all parameters of TPM operation. 

The TpmInput.json file is a basic parameter case with a sample number of 3 and divided into 6 tissue types.
Jsons sample：See TpmInput.json

Folder description:
Eigen-3.4.0: Eigen dependency package 

ITK-prefix: ITK header files 

Sitk: SITK header files 

Libjson: jsonlibjson header files 

SPMSRC: SPM native C code 

Struct2x: Struct2x serialization library

