

using System;
using System.Collections.Generic;
using PylonC.NET;
using System.IO;
using System.Runtime.InteropServices;

namespace SimpleGrab
{
    class SimpleGrab
    {
        //static int captureDuration = (120)+7;
        static int watcherFrequency = 15 * 60; // 15 minutes in seconds
        static int captureFrequency = 5;
        static int numGrabs = 2;
        static float startExposure = 250.0f;
        static float exposureMultiplier = 1.0f;
        
        static void Main(string[] args)
        {    
            PYLON_DEVICE_HANDLE hDev = new PYLON_DEVICE_HANDLE(); /* Handle for the pylon device. */
            initDirectories();
            try
            {
                initDevice(hDev);
                Console.WriteLine("Acquisition Timer init");
                System.Timers.Timer tim_capturer = new System.Timers.Timer();

                tim_capturer.Elapsed += delegate { getFrames(hDev); };
                tim_capturer.Interval = 1000 * captureFrequency;
                tim_capturer.Enabled = true;

                Console.WriteLine("Watchdog Timer init");
                System.Timers.Timer tim_watcher = new System.Timers.Timer();

                tim_watcher.Elapsed += delegate { watcherTimer(tim_capturer); };
                tim_watcher.Interval = 1000 * watcherFrequency;
                tim_watcher.Enabled = true;

                
                Console.WriteLine("Watchdog Timer Start at :" + DateTime.Now.ToString());
                tim_watcher.Start();

                while (tim_watcher.Enabled)
                {
                    System.Threading.Thread.Sleep(5000);
                }
                tim_watcher.Enabled = false;
                tim_watcher.Dispose();
                System.Threading.Thread.Sleep(10000);
                Console.WriteLine("Watchdog Timer Stops at :" + DateTime.Now.ToString());

                finalizeDirectories();

                /* Clean up. Close and release the pylon device. */
                Pylon.DeviceClose(hDev);
                Pylon.DestroyDevice(hDev);

                /* Free memory for grabbing. */
                //imgBuf = null;

                Console.Error.WriteLine("\nPress enter to exit.");
                Console.ReadLine();

                /* Shut down the pylon runtime system. Don't call any pylon method after 
                   calling Pylon.Terminate(). */
                Pylon.Terminate();
            }
            catch (Exception e)
            {
                /* Retrieve the error message. */
                string msg = GenApi.GetLastErrorMessage() + "\n" + GenApi.GetLastErrorDetail();
                Console.Error.WriteLine("Exception caught:");
                Console.Error.WriteLine(e.Message);
                if (msg != "\n")
                {
                    Console.Error.WriteLine("Last error message:");
                    Console.Error.WriteLine(msg);
                }

                try
                {
                    if (hDev.IsValid)
                    {
                        /* ... Close and release the pylon device. */
                        if (Pylon.DeviceIsOpen(hDev))
                        {
                            Pylon.DeviceClose(hDev);
                        }
                        Pylon.DestroyDevice(hDev);
                    }
                }
                catch (Exception)
                {
                    /*No further handling here.*/
                }

                Pylon.Terminate();  /* Releases all pylon resources. */

                Console.Error.WriteLine("\nPress enter to exit.");
                Console.ReadLine();

                Environment.Exit(1);
            }
        }

        static void watcherTimer(System.Timers.Timer tim_capturer) {
            DateTime systemTime = DateTime.Now;
            DateTime startTime = DateTime.Today.AddHours(6D);
            startTime.AddMinutes(30D);
            DateTime endTime = DateTime.Today.AddHours(21D);
            endTime.AddMinutes(30D);

            if ((systemTime > startTime) && (systemTime < endTime) && !tim_capturer.Enabled)
            { // we should start the capture in that case...
                Console.WriteLine("Acquisition Timer Start at :" + DateTime.Now.ToString());
                tim_capturer.Start();
            }
            else if ((systemTime < startTime) && (systemTime > endTime) && tim_capturer.Enabled) 
            { // we should terminate the capture in that case
                tim_capturer.Enabled = false;
                tim_capturer.Dispose();
                Console.WriteLine("Acquisition Timer Stops at :" + DateTime.Now.ToString()); 
            }
        }
        static void initDirectories() {
            if (!System.IO.Directory.Exists("Images"))
                System.IO.Directory.CreateDirectory("Images");
            if (!System.IO.Directory.Exists("Images\\Buffer"))
                System.IO.Directory.CreateDirectory("Images\\Buffer");
            bool file_lock = true;

            if (System.IO.Directory.Exists("Images\\Buffer"))
            {
                String[] bufferedFiles = System.IO.Directory.GetFiles("Images\\Buffer");
                if (bufferedFiles.Length != 0)
                {
                    foreach (String tf in bufferedFiles)
                    {
                        String[] ff = tf.Split('\\');
                        String[] dd = ff[2].Split('_');
                        Console.WriteLine("Images\\Buffer\\" + ff[2]);
                        Console.WriteLine("Images\\" + dd[0] + "\\" + dd[1] + "\\" + dd[2] + "\\" + ff[2]);
                        while (true)
                        {
                            try
                            {
                                System.IO.File.Move("Images\\Buffer\\" + ff[2], "Images\\" + dd[0] + "\\" + dd[1] + "\\" + dd[2] + "\\" + ff[2]);
                                break;
                            }
                            catch (System.IO.IOException e)
                            {
                                if (!IsFileLocked(e))
                                    throw;
                                Console.WriteLine("File is used by another process, will try again...");
                            }
                        }
                    }
                }
            }
        }
        static void finalizeDirectories() {
            if (System.IO.Directory.Exists("Images\\Buffer"))
            {
                String[] bufferedFiles = System.IO.Directory.GetFiles("Images\\Buffer");
                if (bufferedFiles.Length != 0)
                {
                    foreach (String tf in bufferedFiles)
                    {
                        String[] ff = tf.Split('\\');
                        String[] dd = ff[2].Split('_');
                        Console.WriteLine("Images\\Buffer\\" + ff[2]);
                        Console.WriteLine("Images\\" + dd[0] + "\\" + dd[1] + "\\" + dd[2] + "\\" + ff[2]);

                        while (true)
                        {
                            try
                            {
                                System.IO.File.Move("Images\\Buffer\\" + ff[2], "Images\\" + dd[0] + "\\" + dd[1] + "\\" + dd[2] + "\\" + ff[2]);
                                break;
                            }
                            catch (System.IO.IOException e)
                            {
                                if (!IsFileLocked(e))
                                    throw;
                                Console.WriteLine("File is used by another process, will try again...");
                            }
                        }
                    }
                }
                System.IO.Directory.Delete("Images\\Buffer");
            }
        }
        static void initDevice(PYLON_DEVICE_HANDLE hDev)
        {
            uint numDevices;    /* Number of available devices. */

            //PylonBuffer<Byte>[] imgBuf = new PylonBuffer<Byte>[numGrabs];  /* Buffer used for grabbing. */
            bool isAvail;
            int i;
#if DEBUG
            /* This is a special debug setting needed only for GigE cameras.
                See 'Building Applications with pylon' in the Programmer's Guide. */
            Environment.SetEnvironmentVariable("PYLON_GIGE_HEARTBEAT", "300000" /*ms*/);
#endif
            /* Before using any pylon methods, the pylon runtime must be initialized. */
            Pylon.Initialize();

            /* Enumerate all camera devices. You must call 
            PylonEnumerateDevices() before creating a device. */
            numDevices = Pylon.EnumerateDevices();

            if (0 == numDevices)
            {
                throw new Exception("No devices found.");
            }

            /* Get a handle for the first device found.  */
            hDev = Pylon.CreateDeviceByIndex(0);

            /* Before using the device, it must be opened. Open it for configuring
            parameters and for grabbing images. */
            Pylon.DeviceOpen(hDev, Pylon.cPylonAccessModeControl | Pylon.cPylonAccessModeStream);

            /* Set the pixel format to Mono8, where gray values will be output as 8 bit values for each pixel. */
            /* ... Check first to see if the device supports the Mono8 format. */
            isAvail = Pylon.DeviceFeatureIsAvailable(hDev, "EnumEntry_PixelFormat_YUV422Packed");

            if (!isAvail)
            {
                /* Feature is not available. */
                throw new Exception("Device doesn't support the selected pixel format.");
            }

            /* ... Set the pixel format to Mono8. */
            Pylon.DeviceFeatureFromString(hDev, "PixelFormat", "YUV422Packed");


            isAvail = Pylon.DeviceFeatureIsWritable(hDev, "Width") && Pylon.DeviceFeatureIsWritable(hDev, "Height") && Pylon.DeviceFeatureIsWritable(hDev, "OffsetX") && Pylon.DeviceFeatureIsWritable(hDev, "OffsetY");
            // the perimeter of the active pixel area is established here...
            if (isAvail)
            {
                //Pylon.DeviceSetIntegerFeature(hDev, "Width", 1666);
                //Pylon.DeviceSetIntegerFeature(hDev, "Height", 1666);
                //Pylon.DeviceSetIntegerFeature(hDev, "OffsetX", 190);
                //Pylon.DeviceSetIntegerFeature(hDev, "OffsetY", 198);
                Pylon.DeviceSetIntegerFeature(hDev, "Width", 1800);
                Pylon.DeviceSetIntegerFeature(hDev, "Height", 1800);
                Pylon.DeviceSetIntegerFeature(hDev, "OffsetX", 130);
                Pylon.DeviceSetIntegerFeature(hDev, "OffsetY", 130);
            }

            // setings for the capture color balance... 
            isAvail = Pylon.DeviceFeatureIsAvailable(hDev, "EnumEntry_ExposureMode_Timed");
            if (isAvail)
            {
                Pylon.DeviceFeatureFromString(hDev, "ExposureMode", "Timed");
            }

            isAvail = Pylon.DeviceFeatureIsAvailable(hDev, "EnumEntry_ExposureAuto_Off");
            if (isAvail)
            {
                Pylon.DeviceFeatureFromString(hDev, "ExposureAuto", "Off");
            }

            isAvail = Pylon.DeviceFeatureIsWritable(hDev, "ExposureTimeAbs");
            if (isAvail)
            {
                Pylon.DeviceSetFloatFeature(hDev, "ExposureTimeAbs", 300.0);
            }

            // ISO sensitivity of the sensor is arranged in this part
            isAvail = Pylon.DeviceFeatureIsAvailable(hDev, "EnumEntry_GainAuto_Off");
            if (isAvail)
            {
                Pylon.DeviceFeatureFromString(hDev, "GainAuto", "Off");
            }

            isAvail = Pylon.DeviceFeatureIsAvailable(hDev, "EnumEntry_GainSelector_All");
            if (isAvail)
            {
                Pylon.DeviceFeatureFromString(hDev, "GainSelector", "All");
            }
            // here we set the ISO value to be used for the current setting
            isAvail = Pylon.DeviceFeatureIsWritable(hDev, "GainRaw");
            if (isAvail)
            {
                Pylon.DeviceSetIntegerFeature(hDev, "GainRaw", 36);
            }

            // should explore these things here...
            /* Disable acquisition start trigger if available. */
            isAvail = Pylon.DeviceFeatureIsAvailable(hDev, "EnumEntry_TriggerSelector_AcquisitionStart");
            if (isAvail)
            {
                Pylon.DeviceFeatureFromString(hDev, "TriggerSelector", "AcquisitionStart");
                Pylon.DeviceFeatureFromString(hDev, "TriggerMode", "Off");
            }

            /* Disable frame start trigger if available */
            isAvail = Pylon.DeviceFeatureIsAvailable(hDev, "EnumEntry_TriggerSelector_FrameStart");
            if (isAvail)
            {
                Pylon.DeviceFeatureFromString(hDev, "TriggerSelector", "FrameStart");
                Pylon.DeviceFeatureFromString(hDev, "TriggerMode", "Off");
            }


            /* For GigE cameras, we recommend increasing the packet size for better 
               performance. If the network adapter supports jumbo frames, set the packet 
               size to a value > 1500, e.g., to 8192. In this sample, we only set the packet size
               to 1500. */
            /* ... Check first to see if the GigE camera packet size parameter is supported 
                and if it is writable. */
            isAvail = Pylon.DeviceFeatureIsWritable(hDev, "GevSCPSPacketSize");

            if (isAvail)
            {
                /* ... The device supports the packet size feature. Set a value. */
                Pylon.DeviceSetIntegerFeature(hDev, "GevSCPSPacketSize", 1000);
            }
        }
        private static bool IsFileLocked(IOException exception)
        {
            int errorCode = Marshal.GetHRForException(exception) & ((1 << 16) - 1);
            return errorCode == 32 || errorCode == 33;
        }
        /* Simple "image processing" function returning the minimum and maximum gray 
        value of an 8 bit gray value image. */
        static void getProcessed(ref Byte[] resultIm, Byte[] imageBuffer, int numComponents)
        {
            if (resultIm == null)
            {
                resultIm = new Byte[imageBuffer.LongLength];
                resultIm.Initialize();
            }

            for (long i = 0; i < resultIm.LongLength; ++i)
            {
                resultIm[i] += (Byte)(imageBuffer[i] / numComponents);
            }
        }

        static void getFrames(PYLON_DEVICE_HANDLE hDev)
        {
            //const int numGrabs = 4; /* Number of images to grab. */
            /* Grab some images in a loop. */
            Boolean isAvail;

            PylonBuffer<Byte> imgBuf = null;  /* Pylon buffer used for grabbing. */
            PylonBuffer<Byte> resBuf = null;  /* Pylon buffer used for grabbing. */
            Byte[] resByteBuf = null;  /* Buffer used for conversion. */

            PylonGrabResult_t grabResult = null;



            int i;
            for (i = 0; i < numGrabs; ++i)
            {
                isAvail = Pylon.DeviceFeatureIsWritable(hDev, "ExposureTimeAbs");
                if (isAvail)
                {
                    Pylon.DeviceSetFloatFeature(hDev, "ExposureTimeAbs", exposureMultiplier * (i + 1) * startExposure);
                }

                /* Grab one single frame from stream channel 0. The 
                camera is set to "single frame" acquisition mode.
                Wait up to 3000 ms for the image to be grabbed. 
                If imgBuf is null a buffer is automatically created with the right size.*/
                if (!Pylon.DeviceGrabSingleFrame(hDev, 0, ref imgBuf, out grabResult, (uint)(captureFrequency*1000/numGrabs)))
                {
                    /* Timeout occurred. */
                    Console.WriteLine("Frame {0}: timeout.", i + 1);
                }

                /* Check to see if the image was grabbed successfully. */
                if (grabResult.Status == EPylonGrabStatus.Grabbed)
                {
                    /* Success. Perform image processing. */
                    getProcessed(ref resByteBuf, imgBuf.Array, numGrabs);
                    //Console.WriteLine("Grabbed frame {0}. Min. gray value = {1}, Max. gray value = {2}", i + 1, min, max);
                    /* Display image */
                    //Pylon.ImageWindowDisplayImage<Byte>(0, imgBuf, grabResult);

                    /* If we want to save all the components seperately as well... */
                    /*String dir1Name = DateTime.Now.Year.ToString();
                    String dir2Name = DateTime.Now.Month.ToString().PadLeft(2, '0');
                    String dir3Name = DateTime.Now.Day.ToString().PadLeft(2, '0');
                    String fileName = dir1Name + "_" + dir2Name + "_" + dir3Name + "_" + DateTime.Now.Hour.ToString().PadLeft(2, '0') + "_" + DateTime.Now.Minute.ToString().PadLeft(2, '0') + "_" + DateTime.Now.Second.ToString().PadLeft(2, '0') + "_" + DateTime.Now.Millisecond.ToString();


                    if (!System.IO.Directory.Exists("Images\\" + dir1Name))
                        System.IO.Directory.CreateDirectory("Images\\" + dir1Name);
                    if (!System.IO.Directory.Exists("Images\\" + dir1Name + "\\" + dir2Name))
                        System.IO.Directory.CreateDirectory("Images\\" + dir1Name + "\\" + dir2Name);
                    if (!System.IO.Directory.Exists("Images\\" + dir1Name + "\\" + dir2Name + "\\" + dir3Name))
                        System.IO.Directory.CreateDirectory("Images\\" + dir1Name + "\\" + dir2Name + "\\" + dir3Name);
                    Pylon.ImagePersistenceSave<Byte>(EPylonImageFileFormat.ImageFileFormat_Jpeg, "Images\\Buffer\\" + fileName + ".jpeg", imgBuf, EPylonPixelType.PixelType_YUV422packed, (uint)grabResult.SizeX, (uint)grabResult.SizeY, (uint)grabResult.PaddingX, EPylonImageOrientation.ImageOrientation_TopDown);*/
                }
                else if (grabResult.Status == EPylonGrabStatus.Failed)
                {
                    Console.Error.WriteLine("Frame {0} wasn't grabbed successfully.  Error code = {1}", i + 1, grabResult.ErrorCode);
                }
                //System.Threading.Thread.Sleep(500);
            }

            if (System.IO.Directory.Exists("Images\\Buffer"))
            {
                String[] bufferedFiles = System.IO.Directory.GetFiles("Images\\Buffer");
                if (bufferedFiles.Length != 0)
                {
                    foreach (String tf in bufferedFiles)
                    {
                        String[] ff = tf.Split('\\');
                        String[] dd = ff[2].Split('_');
                        Console.WriteLine("Images\\Buffer\\" + ff[2]);
                        Console.WriteLine("Images\\" + dd[0] + "\\" + dd[1] + "\\" + dd[2] + "\\" + ff[2]);
                        while (true)
                        {
                            try
                            {
                                System.IO.File.Move("Images\\Buffer\\" + ff[2], "Images\\" + dd[0] + "\\" + dd[1] + "\\" + dd[2] + "\\" + ff[2]);
                                break;
                            }
                            catch (System.IO.IOException e)
                            {
                                if (!IsFileLocked(e))
                                    throw;
                                Console.WriteLine("File is used by another process, will try again...");
                            }
                        }
                        //System.IO.File.Move("Images\\Buffer\\" + ff[2], "Images\\" + dd[0] + "\\" + dd[1] + "\\" + dd[2] + "\\" + ff[2]);
                    }
                }
            }

            /* Check to see if all in all the image was grabbed successfully. */
            if (grabResult.Status == EPylonGrabStatus.Grabbed)
            {
                resBuf = new PylonBuffer<byte>(resByteBuf);
                //Pylon.ImageWindowDisplayImage<Byte>((uint)i, resBuf, grabResult);
                String dir1Name = DateTime.Now.Year.ToString();
                String dir2Name = DateTime.Now.Month.ToString().PadLeft(2, '0');
                String dir3Name = DateTime.Now.Day.ToString().PadLeft(2, '0');
                String fileName = dir1Name + "_" + dir2Name + "_" + dir3Name + "_" + DateTime.Now.Hour.ToString().PadLeft(2, '0') + "_" + DateTime.Now.Minute.ToString().PadLeft(2, '0') + "_" + DateTime.Now.Second.ToString().PadLeft(2, '0') + "_" + DateTime.Now.Millisecond.ToString();


                if (!System.IO.Directory.Exists("Images\\" + dir1Name))
                    System.IO.Directory.CreateDirectory("Images\\" + dir1Name);
                if (!System.IO.Directory.Exists("Images\\" + dir1Name + "\\" + dir2Name))
                    System.IO.Directory.CreateDirectory("Images\\" + dir1Name + "\\" + dir2Name);
                if (!System.IO.Directory.Exists("Images\\" + dir1Name + "\\" + dir2Name + "\\" + dir3Name))
                    System.IO.Directory.CreateDirectory("Images\\" + dir1Name + "\\" + dir2Name + "\\" + dir3Name);

                Pylon.ImagePersistenceSave<Byte>(EPylonImageFileFormat.ImageFileFormat_Jpeg, "Images\\Buffer\\" + fileName + ".jpeg", imgBuf, EPylonPixelType.PixelType_YUV422packed, (uint)grabResult.SizeX, (uint)grabResult.SizeY, (uint)grabResult.PaddingX, EPylonImageOrientation.ImageOrientation_TopDown);
                /* Display image */
                //Pylon.ImageWindowDisplayImage<Byte>(0, imgBuf, grabResult);
            }



            Console.WriteLine("Grab Procedure Terminates at: " + DateTime.Now.ToString());
        }
    }
}
