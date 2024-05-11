# TPM++
TPM++ Read Me

TPM++ is an advanced C++ implementation of the Tissue Probability Map (TPM) originally found in Statistical Parametric Mapping(SPM). It leverages the Eigen library for efficient matrix operations, replacing the Matlab framework, and incorporates SimpleITK for image processing. The image format and coordinate space are fully compatible with the ITK standards. 

Licensed under GPL v2, TPM++ not only replicates the results of SPM12 but also significantly enhances computational performance. It is designed to operate N times faster than SPM12, where N represents the number of physical CPU cores. 

In its default configuration, TPM++ opts not to utilize hyper-threading. This decision is based on the observation that hyper-threading no longer contributes to substantial performance improvements.

TPM++ also includes a command-line program for calculating Wrap and executing image segmentation, providing a versatile toolset for users.

## Usage

TPM++ provides a command-line program to calculate Wrap and perform image segmentation. The command is as follows:

```bash
TPM++.exe -in image.nii -tpm D:/TPM.nii -out D:/_out.nii -json D:/tpminput.json -AffineOnly 0
```

## Test Case

Input:

![Alt Text](https://github.com/datouxyz/TPMPlusPlus/blob/master/ReadMe.files/ReadMe3193.png)

Outputs:

![Alt Text](https://github.com/datouxyz/TPMPlusPlus/blob/master/ReadMe.files/ReadMe3251.png)
![Alt Text](https://github.com/datouxyz/TPMPlusPlus/blob/master/ReadMe.files/ReadMe3253.png)
![Alt Text](https://github.com/datouxyz/TPMPlusPlus/blob/master/ReadMe.files/ReadMe3274.png)
![Alt Text](https://github.com/datouxyz/TPMPlusPlus/blob/master/ReadMe.files/ReadMe3276.png)

## Parameters

- `-in`: Input image path
- `-tpm`: TPM template path
- `-out`: Output file path, the actual output file will add a suffix to this path:
  - `_FieldInverse`: Inverse deformation field
  - `_FieldForward`: Forward deformation field
  - `_BiasCorrected`: After bias correction
  - `_BiasField`: Bias field
  - `_TissNative`: Segmentation result (in original space)
  - `_TissImported`: Segmentation result (in Dartel space)
  - `_TissModulate`: Segmentation result (in imported space, after modulation)
  - `_TissUnModulate`: Segmentation result (in imported space, before modulation)
- `-AffineOnly`: When set to 1, only the affine transformation is calculated and the affine matrix is output.
- `-json`: The json file that records all parameters of TPM operation.

The TpmInput.json file is a basic parameter case with a sample number of 3 and divided into 6 tissue types.

## JSON Sample

```json
{    
    "m_tissues": [ 
        {
            "m_ngaus": 1,
            "m_NativeTissue_NativeSpace": true,
            "m_NativeTissue_DartelImported": true,
            "m_WarpTissue_Modulated": true,
            "m_WarpTissue_UnModulated": true
        },
        {
            "m_ngaus": 1,
            "m_NativeTissue_NativeSpace": true,
            "m_NativeTissue_DartelImported": true,
            "m_WarpTissue_Modulated": true,
            "m_WarpTissue_UnModulated": true
        },
        {
            "m_ngaus": 2,
            "m_NativeTissue_NativeSpace": true,
            "m_NativeTissue_DartelImported": true,
            "m_WarpTissue_Modulated": true,
            "m_WarpTissue_UnModulated": true
        },
        {
            "m_ngaus": 3,
            "m_NativeTissue_NativeSpace": true,
            "m_NativeTissue_DartelImported": true,
            "m_WarpTissue_Modulated": true,
            "m_WarpTissue_UnModulated": true
        },
        {
            "m_ngaus": 4,
            "m_NativeTissue_NativeSpace": true,
            "m_NativeTissue_DartelImported": true,
            "m_WarpTissue_Modulated": true,
            "m_WarpTissue_UnModulated": true
        },
        {
            "m_ngaus": 2,
            "m_NativeTissue_NativeSpace": true,
            "m_NativeTissue_DartelImported": true,
            "m_WarpTissue_Modulated": true,
            "m_WarpTissue_UnModulated": true
        }
    ],    
    "m_Channels": [
        {
            "m_Bisareg": 0.001000,
            "m_Biasfwhm": 60.000000,
            "m_biascorrect": true,
            "m_biasfield": true
        }
    ],
    "m_nsample": 3,
    "m_fwhm": 0.000000,
    "m_reg": [
        0.000000,
        0.001000,
        0.500000,
        0.050000,
        0.200000
    ],
    "m_lkp": [
        1.000000,
        2.000000,
        3.000000,
        3.000000,
        4.000000,
        4.000000,
        4.000000,
        5.000000,
        5.000000,
        5.000000,
        5.000000,
        6.000000,
        6.000000
    ],
    "m_ptype": 0,
    "m_dFieldInverse": true,
    "m_dFieldForward": true,
    "m_Mrf": false,
    "m_Cleanup": 1
}
```
Folder description:<br>
Eigen-3.4.0: Eigen dependency package <br>
ITK-prefix: ITK header files <br>
Sitk: SITK header files <br>
Libjson: jsonlibjson header files <br>
SPMSRC: SPM native C code <br>
Struct2x: Struct2x serialization library



