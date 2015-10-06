using System;
using System.Collections.Generic;
using PylonC.NETSupportLibrary;
using PylonC.NET;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using ImageAcquisitionLib;
using System.Threading;
using System.IO;

namespace ImageAcquisitionTool_GUI
{
    public partial class ImageAcquisitionGUI : Form
    {
        private bool _isActive;
        private PYLON_IMAGE_FORMAT_CONVERTER_HANDLE _hConverter; //The format converter is used mainly for coverting color images. It is not used for Mono8 or RGBA8packed images.
        private bool _converterOutputFormatIsColor; //The output format of the format converter.
        private PylonGrabResult_t _grabResult;
        private Bitmap _outImage;
        private PylonBuffer<Byte> _pylnBuffer;
        //private Thread _thCapturer;
        ImageAcquisitionTool _ll;
        public ImageAcquisitionGUI()
        {
            InitializeComponent();
            _isActive = true;
            _hConverter = new PYLON_IMAGE_FORMAT_CONVERTER_HANDLE();
            _converterOutputFormatIsColor = false;
            if (_isActive)
            {
                _ll = new ImageAcquisitionTool();
                _ll.statusChangeEvent += new EventHandler(this.statusString_ValueChanged);
                _ll.imageAcquiredEvent += new imageAcquiredEventHandler(this.imageAcquired_ValueChanged);
                _ll.imageGrabbedEvent += new imageGrabbedEventHandler(this.imageGrabbed_ValueChanged);
            }
        }

        private void ImageAcquisitionGUI_Load(object sender, EventArgs e)
        {
            if (_isActive)
            {
                long[] gainLims = _ll.getGainLimits();
                double[] exposureLims = _ll.getExposureLimits();
                sbGain.Maximum = -(int)gainLims[0];
                nmGain.Maximum = (int)gainLims[1];
                sbGain.Minimum = -(int)gainLims[1];
                nmGain.Minimum = (int)gainLims[0];
                sbExp.Maximum = -(int)exposureLims[0];
                nmExp.Maximum = (int)exposureLims[1];
                sbExp.Minimum = -(int)exposureLims[1];
                nmExp.Minimum = (int)exposureLims[0];
            }
            //_thCapturer = new Thread(new ThreadStart(_ll.ProcedureStart));
            pnCamControls1.Enabled = true;
        }

        private void label1_Click(object sender, EventArgs e)
        {

        }

        private void lbExpSt_Click(object sender, EventArgs e)
        {

        }

        private void lbIncGain_Click(object sender, EventArgs e)
        {

        }

        private void btStart_Click(object sender, EventArgs e)
        {
            //_thCapturer.Start();
            if (_isActive)
            {
                _ll.ProcedureStart();
            }
        }

        private void btStop_Click(object sender, EventArgs e)
        {
            if (_isActive)
            {
                _ll.ProcedureStop();
            }
        }

        private void ImageAcquisitionGUI_FormClosed(object sender, FormClosedEventArgs e)
        {
            //_ll.ProcedureTerminate();
        }

        private void ImageAcquisitionGUI_FormClosing(object sender, FormClosingEventArgs e)
        {
            if (_isActive)
            {
                _ll.ProcedureTerminate();
            }
        }

        private void sbGain_Scroll(object sender, ScrollEventArgs e)
        {
            nmGain.Value = -sbGain.Value;
            if (_isActive)
            {
                _ll.setGain((int)nmGain.Value);
            }
        }

        private void nmGain_ValueChanged(object sender, EventArgs e)
        {
            sbGain.Value = -(int)nmGain.Value;
            if (_isActive)
            {
                _ll.setGain((int)nmGain.Value);
            }
        }

        private void sbExp_Scroll(object sender, ScrollEventArgs e)
        {
            nmExp.Value = -sbExp.Value;
            if (_isActive)
            {
                _ll.setStartExposure((int)nmExp.Value);
            }
        }

        private void nmExp_ValueChanged(object sender, EventArgs e)
        {
            sbExp.Value = -(int)nmExp.Value;
            if (_isActive)
            {
                _ll.setStartExposure((int)nmExp.Value);
            }
        }

        private void sbExpInc_Scroll(object sender, ScrollEventArgs e)
        {
            nmExpInc.Value = -sbExpInc.Value;
            if (_isActive)
            {
                _ll.setExposureMultiplier((int)nmExpInc.Value);
            }
        }

        private void nmExpInc_ValueChanged(object sender, EventArgs e)
        {
            sbExpInc.Value = -(int)nmExpInc.Value;
            if (_isActive)
            {
                _ll.setExposureMultiplier((int)nmExpInc.Value);
            }
        }

        private void sbNofShots_Scroll(object sender, ScrollEventArgs e)
        {
            nmNofShots.Value = -sbNofShots.Value;
            if (_isActive)
            {
                _ll.setNumGrabs((int)nmNofShots.Value);
            }
        }

        private void nmNofShots_ValueChanged(object sender, EventArgs e)
        {
            sbNofShots.Value = -(int)nmNofShots.Value;
            if (_isActive)
            {
                _ll.setNumGrabs((int)nmNofShots.Value);
            }
        }

        private void statusString_ValueChanged(object sender, EventArgs e)
        {
            tsStatusLabel.Text = "Status: "+_ll.statusMessage;
        }

        private void imageAcquired_ValueChanged(object sender, imageAcquiredEventArgs e)
        {
            pbImage.Image = e.getImgBuffer();
        }

        private void imageGrabbed_ValueChanged(object sender, imageGrabbedEventArgs e)
        {
            // Reference to the buffer attached to the grab result handle.
            _pylnBuffer = e.getImgBuffer();
            _grabResult = e.getGrabResult();
            PylonBuffer<Byte> convertedBuffer = null;
            
            if (!_hConverter.IsValid)
            {
                _hConverter = Pylon.ImageFormatConverterCreate(); /* Create the converter. */
                _converterOutputFormatIsColor = !Pylon.IsMono(_grabResult.PixelType) || Pylon.IsBayer(_grabResult.PixelType);
            }

            // Perform the conversion. If the buffer is null a new one is automatically created.
            Pylon.ImageFormatConverterSetOutputPixelFormat(_hConverter, _converterOutputFormatIsColor ? EPylonPixelType.PixelType_BGRA8packed : EPylonPixelType.PixelType_Mono8);
            Pylon.ImageFormatConverterConvert(_hConverter, ref convertedBuffer, _pylnBuffer, _grabResult.PixelType, (uint)_grabResult.SizeX, (uint)_grabResult.SizeY, (uint)_grabResult.PaddingX, EPylonImageOrientation.ImageOrientation_TopDown);

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
            }
            //pbImage.Refresh();
            pbImage.Image = _outImage;
        }
    }
}
