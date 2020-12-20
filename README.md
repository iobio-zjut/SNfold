# SNfold
### a sequential niche multimodal conformation sampling algorithm for protein structure prediction



**Developer:**   
                Yuhao Xia
                College of Information Engineering
                Zhejiang University of Technology, Hangzhou 310023, China
		Email: xiayh@zjut.edu.cn

**Contact:**  
                Guijun Zhang, Prof
                College of Information Engineering
                Zhejiang University of Technology, Hangzhou 310023, China
                Email: zgj@zjut.edu.cn

## 1. INSTALLATION
Please Follow the below steps to install and configure SNfold:

- Download Rosetta3.10 source package from https://www.rosettacommons.org/software/, and extract it to ``"~/"`` directory where SNfold package is located.

- Copy and paste ``"ClassicAbinitio.cc"`` and ``"ClassicAbinitio.hh"`` from ``"src/"`` folder in SNfold package to ``"~/Rosetta/main/source/src/protocols/abinitio/"`` folder in Rosetta.

- Copy and paste ``"MonteCarlo.cc"`` and ``"MonteCarlo.hh"`` from ``"src/"`` folder in SNfold package to ``"~/Rosetta/main/source/src/protocols/moves/"`` folder in Rosetta.

- Compile SNfold source code using the following commands:

```
 $ cd ~/Rosetta/main/source/
 $ ./scons.py AbinitioRelax -j<NumOfJobs> mode=release bin
```

## 2. INPUT
SNfold requires four files to generate models:

	-fasta				: fasta file
	-alldist.txt			: distance map file
	-aat000_03_05.200_v1_3		: fragment library with fragment lenth 3
	-aat000_09_05.200_v1_3		: fragment library with fragment lenth 9

## 3. RUN
Please follow the below steps to run SNfold:

- Go to the ``"example/"`` folder of SNfold.

- Run SNfold with the following command:

```
   $ ./run.sh
```

- Five models are generated in the ``"output_files/"`` folder.

## 4. OUTPUT
Output files of SNfold are stored in the ``"example/output_files/"`` folder, including five predicted models.

## 5. DISCLAIMER
The executable software and the source code of SNfold is distributed free of charge 
as it is to any non-commercial users. The authors hold no liabilities to the performance 
of the program.
