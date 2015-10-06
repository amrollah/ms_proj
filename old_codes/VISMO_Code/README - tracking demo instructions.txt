Instructions:


1.- 	Visual Studio 2010 Installation
	You need a software package to load an image disk. I recomend MagicDisc
	Go to Go to Sofware directory (I:\110_Project\AGI\C1-164-12_01 Contr and Design of Conc Solar Power Pl\5_Results\Burak Results\Internship\Software) and load the image 	en_visual_studio_2010_ultimate_x86_dvd_509116.iso
	Visual Studio 2010 should be installed
2.- 	OpenCV Installation
	Go to Sofware directory (I:\110_Project\AGI\C1-164-12_01 Contr and Design of Conc Solar Power Pl\5_Results\Burak Results\Internship\Software) and install OpenCV-2.4.5
	Put the installation in your C:\ directory and name it opencv-2.4.5 (e.g. C:\opencv-2.4.5)
	Once you have completed the installation go to "Computer" --> right click on properties --> Advance System Settings --> environmental variables --> System Variables, click on the variable Path and edit it. 
	Then add the following path "C:\opencv-2.4.5\opencv\build\x64\vc10\bin" to the path variable.
	OpenCV libraries are now installed.
3.-	Go to Sofware directory (I:\110_Project\AGI\C1-164-12_01 Contr and Design of Conc Solar Power Pl\5_Results\Burak Results\Internship\Software) and unzip mexopencv-master.zip
	You need to compile the files. Make sure that you have installed a compiler
	Type "mex -setup" to detect the list of available comiles. I reccomend using the Windows SDK Compiler
	Compile mexopencv by typing "mexopencv.make('opencv_path', 'c:\your\path\to\opencv')"   (Note that if there is a compiled version done with the same MATLAB version as yours, you just need to add them to the path, otherwise add the source to path and do this step above.)
	You can read the readme file: README.markwdown if in doubt.
	You can execute some of the examples from samples/ to test the installation
4.- 	go to I:\110_Project\AGI\C1-164-12_01 Contr and Design of Conc Solar Power Pl\5_Results\Burak Results\Internship\Packages  and copy the folders "ImageAcquisitionTool_Auto" and "ImageAcquisitionTool_Manual" into your visual studio 2010 project installation (e.g  "C:\Users\your-windows-user-name\Documents\Visual Studio 2010\Projects")
5.-  	go to I:\110_Project\AGI\C1-164-12_01 Contr and Design of Conc Solar Power Pl\5_Results\Burak Results\Internship\Packages and copy the folder "MATLAB" to your matlab working direactory (e.g. "C:\Users\your-windows-user-name\Documents\MATLAB")
6.-	go to I:\110_Project\AGI\C1-164-12_01 Contr and Design of Conc Solar Power Pl\5_Results\Burak Results\Internship\Packages\MATLAB\CH-CR roof images and execute the file "trackingDemo_CamMeas.m" if it works that means you have followed the instructions correctly. 



For questions or help, drop me a line: luis.dominguez@ch.abb.com

	 