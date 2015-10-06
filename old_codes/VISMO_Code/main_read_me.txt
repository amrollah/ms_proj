Hello, 

If you are reading this file right now, I assume that means you are interested in running and understanding my code.
I will happily try to help you with that. 

The entry point of the code is the file named trackingDemo_CamCalib_v2. The parameters to be input to the 
algorithm are in the 'config.xml' file. What each parameter does can be understood following the line where that 
is used to assign the functional variables inside the code in the 'setupSystemObjects' function of the 
'multiObjectTrackingCV_CamCalib' file. 

The most important parameter to be set is the path to the folder that contains images. This path should be set
correctly. Most of the other parameters can assume default values. 

The camera records an image every 3 seconds. The frame rate of the algoritm states how many frames are skipped 
during operation. So if the frame rate is 4 this means there are 12 seconds between two frames input.
The predictions are done for every single frame in the horizon, this means that with every recursive nth prediction
n*12 seconds ahead is predicted. Horizon is the number of these recursive predictions. So with a horizon h, the 
algorithm generates predictions up to h*12 seconds ahead. For instance if the horizon is 25, algorithm predicts
upto 25*12 = 300 seconds (5 minutes) with 12 seconds increments. 

The code also has a second form given in cloudTrackingCV file. It is implemented in an object oriented way, and does
exactly the same thing. This was an effort to ease understanding the code, but one can understand it by reading 
the multiObjectTrackingCV_CamCalib code better. Every line of that code is commented, the design is given in my
thesis. 

Also there is a project called pilot_trial2 in the folder. That is a project for compiling the whole matlab code into
a standalone executable. I did it for trackingDemo_CamCalib and by then the structure was a bit different but essentially
the structure is the same. Just what is needed to be done is, if there is any input data to the algorithm, these should be
put into a single folder, their names should be stated in the config.xml file and "to the algorithm variables, just above 
the imageFolder, another variable called path should be added where it contains the path to all the data".

I hope I told the story well enough (I tried to be as thorough as possible as if I am telling the story to 
Bilal :P), but there are many other supplementary tools of the system such as labelers and calibration tools.
Anything else is required please don't hesitate ask me on burak.zeydan@epfl.ch ...

Good luck and enjoy,

Ciao
Burak Zeydan, 25.07.2014 
ABB Corporate Research Center
Baden Dattwil, Switzerland