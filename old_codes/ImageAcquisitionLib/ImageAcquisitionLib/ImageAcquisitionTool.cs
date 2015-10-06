using System;
using System.Drawing;
using System.Collections.Generic;
using PylonC.NET;
using PylonC.NETSupportLibrary;
using System.IO;
using System.Linq;
using System.Text;
using System.Runtime.InteropServices;

namespace ImageAcquisitionLib
{
    public class ImageAcquisitionTool
    {
        //static int captureDuration = (120)+7;
        public event EventHandler statusChangeEvent;
        public event imageAcquiredEventHandler imageAcquiredEvent;
        public event imageGrabbedEventHandler imageGrabbedEvent;
        private int _watcherFrequency;
        private int _captureFrequency;
        private int _gain;
        private int _numGrabs;
        private float _startExposure;
        private float _exposureMultiplier;
        private double _minExposure;
        private double _maxExposure;
        private long _minGain;
        private long _maxGain;
        private bool _systemEnabled;
        private Object _thisLock;
        private Bitmap _outImage;

        public string statusMessage;
        private System.Timers.Timer _timCapturer;
        private System.Timers.Timer _timWatcher;

        private PylonBuffer<Byte> _imgBuf;  /* Pylon buffer used for grabbing. */
        private PylonBuffer<Byte> _resBuf;  /* Pylon buffer used for grabbing. */
        private Byte[] _resByteBuf;  /* Buffer used for conversion. */
        private PylonGrabResult_t _grabResult;
        private PYLON_IMAGE_FORMAT_CONVERTER_HANDLE _hConverter; /* The format converter is used mainly for coverting color images. It is not used for Mono8 or RGBA8packed images. */
        private bool _converterOutputFormatIsColor;/* The output format of the format converter. */
        private PYLON_DEVICE_HANDLE _hDev; /* Handle for the pylon device. */
        
        public ImageAcquisitionTool() {
            _watcherFrequency = 5;//15 * 60; // 15 minutes in seconds
            _captureFrequency = 3;
            _numGrabs = 2;
            _gain = 50;
            _startExposure = 250.0f;
            _exposureMultiplier = 1.0f;
            _systemEnabled = false;
            _thisLock = new Object();
            statusMessage = string.Empty;

            _imgBuf = null;  /* Pylon buffer used for grabbing. */
            _resBuf = null;  /* Pylon buffer used for grabbing. */
            _resByteBuf = null;  /* Buffer used for conversion. */
            _grabResult = null;
            _converterOutputFormatIsColor = false;
            _outImage = null;

            _hDev = new PYLON_DEVICE_HANDLE(); /* Handle for the pylon device. */
            _hConverter = new PYLON_IMAGE_FORMAT_CONVERTER_HANDLE();
            initDirectories();

            try {
                
                initDevice(out _hDev);
                _minExposure = Pylon.DeviceGetFloatFeatureMin(_hDev, "ExposureTimeAbs");
                _maxExposure = Pylon.DeviceGetFloatFeatureMax(_hDev, "ExposureTimeAbs");
                _minGain = Pylon.DeviceGetIntegerFeatureMin(_hDev, "GainRaw");
                _maxGain = Pylon.DeviceGetIntegerFeatureMax(_hDev, "GainRaw");

                //setStatusMessage("Acquisition Timer init");
                _timCapturer = new System.Timers.Timer();

                _timCapturer.Elapsed += delegate { getFrames(_hDev); };
                _timCapturer.Interval = 1000 * _captureFrequency;
                //_timCapturer.Enabled = false;

                setStatusMessage("Watchdog Timer init");
                _timWatcher = new System.Timers.Timer();

                _timWatcher.Elapsed += delegate { watcherTimer(_timCapturer, ref _systemEnabled); };
                _timWatcher.Interval = 1000 * _watcherFrequency;
                //_timWatcher.Enabled = true;
                setStatusMessage("Watchdog Timer Start at :" + DateTime.Now.ToString());
                _timWatcher.Start();
            }
            catch (Exception e) {
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
                    if (_hDev.IsValid)
                    {
                        /* ... Close and release the pylon device. */
                        if (Pylon.DeviceIsOpen(_hDev))
                        {
                            Pylon.DeviceClose(_hDev);
                        }
                        Pylon.DestroyDevice(_hDev);
                    }
                }
                catch (Exception)
                {
                    /*No further handling here.*/
                }

                Pylon.Terminate();  /* Releases all pylon resources. */

                //Console.Error.WriteLine("\nPress enter to exit.");
                //Console.ReadLine();

                //Environment.Exit(1);
            }
        }

        public void ProcedureStart()
        { 
            lock (_thisLock)
            {
                _systemEnabled = true;
            }
        }

        public void ProcedureStop()
        {
            lock (_thisLock)
            {
                _systemEnabled = false;
            } 
        }

        public void ProcedureTerminate()
        {
            lock(_thisLock)
            {
                _systemEnabled = false;
            }
            System.Threading.Thread.Sleep(5000);
            _timWatcher.Enabled = false;
            _timWatcher.Dispose();
            System.Threading.Thread.Sleep(5000);
            setStatusMessage("Watchdog Timer Terminates at :" + DateTime.Now.ToString());

            finalizeDirectories();

            /* Clean up. Close and release the pylon device. */
            Pylon.DeviceClose(_hDev);
            Pylon.DestroyDevice(_hDev);
            /* Shut down the pylon runtime system. Don't call any pylon method after 
                   calling Pylon.Terminate(). */
            Pylon.Terminate();
        }
        
        public void setWatcherFrequency(int watcherFrequency) {
            _watcherFrequency = watcherFrequency;
        }
        public int getWatcherFrequency() {
            return _watcherFrequency;
        }
        public void setCaptureFrequency(int captureFrequency)
        {
            _captureFrequency = captureFrequency;
        }
        public int getCaptureFrequency()
        {
            return _captureFrequency;
        }
        public void setNumGrabs(int numGrabs)
        {
            _numGrabs = numGrabs;
        }
        public int getNumGrabs()
        {
            return _numGrabs;
        }
        public void setStartExposure(float startExposure)
        {
            _startExposure = startExposure;
        }
        public float getStartExposure()
        {
            return _startExposure;
        }
        public void setExposureMultiplier(float exposureMultiplier)
        {
            _exposureMultiplier = exposureMultiplier;
        }
        public float getExposureMultiplier()
        {
            return _exposureMultiplier;
        }
        public double[] getExposureLimits()
        {
            return new double[] {_minExposure, _maxExposure};
        }
        public void setGain(int gain)
        {
            _gain = gain;
        }
        public int getGain()
        {
            return _gain;
        }
        public long[] getGainLimits() { 
            return new long[] {_minGain, _maxGain};
        }
        public Byte[] getImage() 
        {
            return _resByteBuf;
        }

        private void setStatusMessage(string s) {
            if (statusMessage.Length != 0)
            {
                statusMessage = string.Empty;
            }
                statusMessage += s;
                EventArgs args = new EventArgs();
                OnStatusChangeEvent(args);
        }
        private void watcherTimer(System.Timers.Timer tim_capturer, ref bool systemEnabled)
        {
            DateTime systemTime = DateTime.Now;
            DateTime startTime = DateTime.Today.AddHours(6D);
            startTime.AddMinutes(30D);
            DateTime endTime = DateTime.Today.AddHours(21D);
            endTime.AddMinutes(30D);
            //setStatusMessage("Hey There : " + systemEnabled.ToString());
            if ((systemTime > startTime) && (systemTime < endTime) && !tim_capturer.Enabled && systemEnabled)
            { // we should start the capture in that case...
                setStatusMessage("Acquisition Timer Start at :" + DateTime.Now.ToString());
                tim_capturer.Start();
            }
            else if (((systemTime < startTime) && (systemTime > endTime) && tim_capturer.Enabled) || (!systemEnabled && tim_capturer.Enabled))
            { // we should terminate the capture in that case
                tim_capturer.Enabled = false;
                //tim_capturer.Dispose();
                tim_capturer.Stop();
                setStatusMessage("Acquisition Timer Stops at :" + DateTime.Now.ToString());
            }
        }
        private void initDirectories()
        {
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
                        setStatusMessage("Images\\Buffer\\" + ff[2]);
                        setStatusMessage("Images\\" + dd[0] + "\\" + dd[1] + "\\" + dd[2] + "\\" + ff[2]);
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
                                setStatusMessage("File is used by another process, will try again...");
                            }
                        }
                    }
                }
            }
        }
        private void finalizeDirectories()
        {
            if (System.IO.Directory.Exists("Images\\Buffer"))
            {
                String[] bufferedFiles = System.IO.Directory.GetFiles("Images\\Buffer");
                if (bufferedFiles.Length != 0)
                {
                    foreach (String tf in bufferedFiles)
                    {
                        String[] ff = tf.Split('\\');
                        String[] dd = ff[2].Split('_');
                        setStatusMessage("Images\\Buffer\\" + ff[2]);
                        setStatusMessage("Images\\" + dd[0] + "\\" + dd[1] + "\\" + dd[2] + "\\" + ff[2]);

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
                                setStatusMessage("File is used by another process, will try again...");
                            }
                        }
                    }
                }
                System.IO.Directory.Delete("Images\\Buffer");
            }
        }
        private void initDevice(out PYLON_DEVICE_HANDLE hDev)
        {
            uint numDevices;    /* Number of available devices. */

            //PylonBuffer<Byte>[] imgBuf = new PylonBuffer<Byte>[numGrabs];  /* Buffer used for grabbing. */
            bool isAvail;
            //int i;
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
                Pylon.DeviceSetIntegerFeature(hDev, "GainRaw", _gain);
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
        private bool IsFileLocked(IOException exception)
        {
            int errorCode = Marshal.GetHRForException(exception) & ((1 << 16) - 1);
            return errorCode == 32 || errorCode == 33;
        }
        /* Simple "image processing" function returning the minimum and maximum gray 
        value of an 8 bit gray value image. */
        private void getProcessed(ref Byte[] resultIm, Byte[] imageBuffer, int numComponents)
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
        private void getFrames(PYLON_DEVICE_HANDLE hDev)
        {
            //const int numGrabs = 4; /* Number of images to grab. */
            /* Grab some images in a loop. */
            Boolean isAvail;

            //PylonBuffer<Byte> imgBuf = null;  /* Pylon buffer used for grabbing. */
            //PylonBuffer<Byte> resBuf = null;  /* Pylon buffer used for grabbing. */
            //Byte[] resByteBuf = null;  /* Buffer used for conversion. */

            //PylonGrabResult_t grabResult = null;
            _imgBuf = null;  /* Pylon buffer used for grabbing. */
            _resBuf = null;  /* Pylon buffer used for grabbing. */
            _resByteBuf = null;  /* Buffer used for conversion. */
            lock (_thisLock)
            {
                _grabResult = null;
            }

            isAvail = Pylon.DeviceFeatureIsWritable(_hDev, "GainRaw");
            if (isAvail)
            {
                Pylon.DeviceSetIntegerFeature(hDev, "GainRaw", _gain);
            }

            int i;
            for (i = 0; i < _numGrabs; ++i)
            {
                isAvail = Pylon.DeviceFeatureIsWritable(hDev, "ExposureTimeAbs");
                if (isAvail)
                {
                    Pylon.DeviceSetFloatFeature(hDev, "ExposureTimeAbs", _exposureMultiplier * (i + 1) * _startExposure);
                }

                /* Grab one single frame from stream channel 0. The 
                camera is set to "single frame" acquisition mode.
                Wait up to 3000 ms for the image to be grabbed. 
                If imgBuf is null a buffer is automatically created with the right size.*/
                if (!Pylon.DeviceGrabSingleFrame(hDev, 0, ref _imgBuf, out _grabResult, (uint)(_captureFrequency * 1000 / _numGrabs)))
                {
                    /* Timeout occurred. */
                    setStatusMessage("Frame" + (i+1).ToString() + "timeout.");
                }

                /* Check to see if the image was grabbed successfully. */
                if (_grabResult.Status == EPylonGrabStatus.Grabbed)
                {
                    /* Success. Perform image processing. */
                    getProcessed(ref _resByteBuf, _imgBuf.Array, _numGrabs);
                    //setStatusMessage("Grabbed frame {0}. Min. gray value = {1}, Max. gray value = {2}", i + 1, min, max);
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
                else if (_grabResult.Status == EPylonGrabStatus.Failed)
                {
                    Console.Error.WriteLine("Frame {0} wasn't grabbed successfully.  Error code = {1}", i + 1, _grabResult.ErrorCode);
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
                        setStatusMessage("Images\\Buffer\\" + ff[2]);
                        setStatusMessage("Images\\" + dd[0] + "\\" + dd[1] + "\\" + dd[2] + "\\" + ff[2]);
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
                                setStatusMessage("File is used by another process, will try again...");
                            }
                        }
                        //System.IO.File.Move("Images\\Buffer\\" + ff[2], "Images\\" + dd[0] + "\\" + dd[1] + "\\" + dd[2] + "\\" + ff[2]);
                    }
                }
            }

            /* Check to see if all in all the image was grabbed successfully. */
            if (_grabResult.Status == EPylonGrabStatus.Grabbed)
            {
                _resBuf = new PylonBuffer<byte>(_resByteBuf);
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

                Pylon.ImagePersistenceSave<Byte>(EPylonImageFileFormat.ImageFileFormat_Jpeg, "Images\\Buffer\\" + fileName + ".jpeg", _resBuf, EPylonPixelType.PixelType_YUV422packed, (uint)_grabResult.SizeX, (uint)_grabResult.SizeY, (uint)_grabResult.PaddingX, EPylonImageOrientation.ImageOrientation_TopDown);
                //Display image
                //Pylon.ImageWindowDisplayImage<Byte>(0, imgBuf, grabResult);
                
                // Create a new format converter if needed.
                /*if (!_hConverter.IsValid)
                {
                    _hConverter = Pylon.ImageFormatConverterCreate(); // Create the converter.
                    _converterOutputFormatIsColor = !Pylon.IsMono(_grabResult.PixelType) || Pylon.IsBayer(_grabResult.PixelType);
                }

                // Reference to the buffer attached to the grab result handle.
                PylonBuffer<Byte> convertedBuffer = null;

                // Perform the conversion. If the buffer is null a new one is automatically created.
                Pylon.ImageFormatConverterSetOutputPixelFormat(_hConverter, _converterOutputFormatIsColor ? EPylonPixelType.PixelType_BGRA8packed : EPylonPixelType.PixelType_Mono8);
                Pylon.ImageFormatConverterConvert(_hConverter, ref convertedBuffer, _resBuf, _grabResult.PixelType, (uint)_grabResult.SizeX, (uint)_grabResult.SizeY, (uint)_grabResult.PaddingX, EPylonImageOrientation.ImageOrientation_TopDown);
                
                // Add the image data.
                //_outImage = new Image(_grabResult.SizeX, _grabResult.SizeY, convertedBuffer.Array, _converterOutputFormatIsColor);
                if (BitmapFactory.IsCompatible(_outImage, _grabResult.SizeX, _grabResult.SizeY, _converterOutputFormatIsColor))
                {
                    BitmapFactory.UpdateBitmap(_outImage, convertedBuffer.Array, _grabResult.SizeX, _grabResult.SizeY, _converterOutputFormatIsColor);
                }
                else 
                {
                    BitmapFactory.CreateBitmap(out _outImage, _grabResult.SizeX, _grabResult.SizeY, _converterOutputFormatIsColor);
                    BitmapFactory.UpdateBitmap(_outImage, convertedBuffer.Array, _grabResult.SizeX, _grabResult.SizeY, _converterOutputFormatIsColor);
                }*/
                //imageAcquiredEventArgs args = new imageAcquiredEventArgs(_outImage);
                //OnImageAcquiredEvent(args);
                imageGrabbedEventArgs args = new imageGrabbedEventArgs(_resBuf,_grabResult);
                OnImageGrabbedEvent(args);
            }
            setStatusMessage("Grab Procedure Terminates at: " + DateTime.Now.ToString());
        }
        protected void OnStatusChangeEvent(EventArgs e)
        {
            EventHandler handler = statusChangeEvent;
            if (handler != null)
            {
                handler(this, e);
            }
        }
        protected void OnImageAcquiredEvent(imageAcquiredEventArgs e)
        {
            imageAcquiredEventHandler handler = imageAcquiredEvent;
            if (handler != null)
            {
                handler(this, e);
            }
        }

        protected void OnImageGrabbedEvent(imageGrabbedEventArgs e)
        {
            imageGrabbedEventHandler handler = imageGrabbedEvent;
            if (handler != null)
            {
                handler(this, e);
            }
        }
    }

    public delegate void imageAcquiredEventHandler(object sender, imageAcquiredEventArgs e);

    public delegate void imageGrabbedEventHandler(object sender, imageGrabbedEventArgs e);
    
    public class imageAcquiredEventArgs : EventArgs 
    {
        private Bitmap _imgBuffer;

        public imageAcquiredEventArgs(Bitmap imgBuffer)
        {
            //_imgBuffer = new Byte[imgBuffer.Length];
            //Buffer.BlockCopy(imgBuffer, 0, _imgBuffer, 0, imgBuffer.Length);
            _imgBuffer = imgBuffer;
        }

        public Bitmap getImgBuffer() 
        {
            return _imgBuffer;
        }
    }

    public class imageGrabbedEventArgs : EventArgs
    {
        private PylonBuffer<Byte> _imgBuffer;
        private PylonGrabResult_t _grabResult;
        public imageGrabbedEventArgs(PylonBuffer<Byte> imgBuffer, PylonGrabResult_t grabResult)
        {
            _imgBuffer = new PylonBuffer<Byte>(imgBuffer.Array);
            _grabResult = grabResult;
        }

        public PylonBuffer<Byte> getImgBuffer()
        {
            return _imgBuffer;
        }

        public PylonGrabResult_t getGrabResult() 
        {
            return _grabResult;
        }
    }
}
