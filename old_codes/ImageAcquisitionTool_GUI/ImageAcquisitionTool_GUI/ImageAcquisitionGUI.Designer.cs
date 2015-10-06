namespace ImageAcquisitionTool_GUI
{
    partial class ImageAcquisitionGUI
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.pbImage = new System.Windows.Forms.PictureBox();
            this.sbGain = new System.Windows.Forms.VScrollBar();
            this.nmGain = new System.Windows.Forms.NumericUpDown();
            this.nmExp = new System.Windows.Forms.NumericUpDown();
            this.sbExp = new System.Windows.Forms.VScrollBar();
            this.nmExpLow = new System.Windows.Forms.NumericUpDown();
            this.sbExpLow = new System.Windows.Forms.VScrollBar();
            this.nmExpHigh = new System.Windows.Forms.NumericUpDown();
            this.sbExpHigh = new System.Windows.Forms.VScrollBar();
            this.lbGain = new System.Windows.Forms.Label();
            this.lbExp = new System.Windows.Forms.Label();
            this.lbExpLow = new System.Windows.Forms.Label();
            this.lbExpHigh = new System.Windows.Forms.Label();
            this.pnCamControls1 = new System.Windows.Forms.Panel();
            this.lbNofShots = new System.Windows.Forms.Label();
            this.nmNofShots = new System.Windows.Forms.NumericUpDown();
            this.sbNofShots = new System.Windows.Forms.VScrollBar();
            this.lbExpInc = new System.Windows.Forms.Label();
            this.nmExpInc = new System.Windows.Forms.NumericUpDown();
            this.sbExpInc = new System.Windows.Forms.VScrollBar();
            this.btStart = new System.Windows.Forms.Button();
            this.btStop = new System.Windows.Forms.Button();
            this.statusStrip1 = new System.Windows.Forms.StatusStrip();
            this.tsStatusLabel = new System.Windows.Forms.ToolStripStatusLabel();
            ((System.ComponentModel.ISupportInitialize)(this.pbImage)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nmGain)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nmExp)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nmExpLow)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nmExpHigh)).BeginInit();
            this.pnCamControls1.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.nmNofShots)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.nmExpInc)).BeginInit();
            this.statusStrip1.SuspendLayout();
            this.SuspendLayout();
            // 
            // pbImage
            // 
            this.pbImage.Location = new System.Drawing.Point(451, 3);
            this.pbImage.Name = "pbImage";
            this.pbImage.Size = new System.Drawing.Size(447, 447);
            this.pbImage.SizeMode = System.Windows.Forms.PictureBoxSizeMode.StretchImage;
            this.pbImage.TabIndex = 0;
            this.pbImage.TabStop = false;
            // 
            // sbGain
            // 
            this.sbGain.Location = new System.Drawing.Point(7, 40);
            this.sbGain.Maximum = 35;
            this.sbGain.Minimum = 1;
            this.sbGain.Name = "sbGain";
            this.sbGain.Size = new System.Drawing.Size(50, 111);
            this.sbGain.TabIndex = 2;
            this.sbGain.Value = 16;
            this.sbGain.Scroll += new System.Windows.Forms.ScrollEventHandler(this.sbGain_Scroll);
            // 
            // nmGain
            // 
            this.nmGain.Location = new System.Drawing.Point(7, 166);
            this.nmGain.Name = "nmGain";
            this.nmGain.Size = new System.Drawing.Size(50, 20);
            this.nmGain.TabIndex = 3;
            this.nmGain.ValueChanged += new System.EventHandler(this.nmGain_ValueChanged);
            // 
            // nmExp
            // 
            this.nmExp.Location = new System.Drawing.Point(70, 166);
            this.nmExp.Name = "nmExp";
            this.nmExp.Size = new System.Drawing.Size(50, 20);
            this.nmExp.TabIndex = 5;
            this.nmExp.ValueChanged += new System.EventHandler(this.nmExp_ValueChanged);
            // 
            // sbExp
            // 
            this.sbExp.Location = new System.Drawing.Point(70, 40);
            this.sbExp.Maximum = 35;
            this.sbExp.Minimum = 1;
            this.sbExp.Name = "sbExp";
            this.sbExp.Size = new System.Drawing.Size(50, 111);
            this.sbExp.TabIndex = 4;
            this.sbExp.Value = 16;
            this.sbExp.Scroll += new System.Windows.Forms.ScrollEventHandler(this.sbExp_Scroll);
            // 
            // nmExpLow
            // 
            this.nmExpLow.Location = new System.Drawing.Point(133, 166);
            this.nmExpLow.Name = "nmExpLow";
            this.nmExpLow.Size = new System.Drawing.Size(50, 20);
            this.nmExpLow.TabIndex = 7;
            // 
            // sbExpLow
            // 
            this.sbExpLow.Location = new System.Drawing.Point(133, 40);
            this.sbExpLow.Maximum = 35;
            this.sbExpLow.Minimum = 1;
            this.sbExpLow.Name = "sbExpLow";
            this.sbExpLow.Size = new System.Drawing.Size(50, 111);
            this.sbExpLow.TabIndex = 6;
            this.sbExpLow.Value = 16;
            // 
            // nmExpHigh
            // 
            this.nmExpHigh.Location = new System.Drawing.Point(196, 166);
            this.nmExpHigh.Name = "nmExpHigh";
            this.nmExpHigh.Size = new System.Drawing.Size(50, 20);
            this.nmExpHigh.TabIndex = 9;
            // 
            // sbExpHigh
            // 
            this.sbExpHigh.Location = new System.Drawing.Point(196, 40);
            this.sbExpHigh.Maximum = 35;
            this.sbExpHigh.Minimum = 1;
            this.sbExpHigh.Name = "sbExpHigh";
            this.sbExpHigh.Size = new System.Drawing.Size(50, 111);
            this.sbExpHigh.TabIndex = 8;
            this.sbExpHigh.Value = 16;
            // 
            // lbGain
            // 
            this.lbGain.AutoSize = true;
            this.lbGain.Font = new System.Drawing.Font("Microsoft Sans Serif", 8F, System.Drawing.FontStyle.Bold, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.lbGain.Location = new System.Drawing.Point(16, 12);
            this.lbGain.Name = "lbGain";
            this.lbGain.Size = new System.Drawing.Size(33, 13);
            this.lbGain.TabIndex = 10;
            this.lbGain.Text = "Gain";
            this.lbGain.TextAlign = System.Drawing.ContentAlignment.MiddleCenter;
            this.lbGain.Click += new System.EventHandler(this.label1_Click);
            // 
            // lbExp
            // 
            this.lbExp.AutoSize = true;
            this.lbExp.Font = new System.Drawing.Font("Microsoft Sans Serif", 8F, System.Drawing.FontStyle.Bold, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.lbExp.Location = new System.Drawing.Point(66, 12);
            this.lbExp.Name = "lbExp";
            this.lbExp.Size = new System.Drawing.Size(59, 13);
            this.lbExp.TabIndex = 11;
            this.lbExp.Text = "Exposure";
            this.lbExp.TextAlign = System.Drawing.ContentAlignment.MiddleCenter;
            this.lbExp.Click += new System.EventHandler(this.lbExpSt_Click);
            // 
            // lbExpLow
            // 
            this.lbExpLow.AutoSize = true;
            this.lbExpLow.Font = new System.Drawing.Font("Microsoft Sans Serif", 8F, System.Drawing.FontStyle.Bold, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.lbExpLow.Location = new System.Drawing.Point(129, 12);
            this.lbExpLow.Name = "lbExpLow";
            this.lbExpLow.Size = new System.Drawing.Size(59, 13);
            this.lbExpLow.TabIndex = 12;
            this.lbExpLow.Text = "Exp. Low";
            this.lbExpLow.TextAlign = System.Drawing.ContentAlignment.MiddleCenter;
            // 
            // lbExpHigh
            // 
            this.lbExpHigh.AutoSize = true;
            this.lbExpHigh.Font = new System.Drawing.Font("Microsoft Sans Serif", 8F, System.Drawing.FontStyle.Bold, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.lbExpHigh.Location = new System.Drawing.Point(190, 12);
            this.lbExpHigh.Name = "lbExpHigh";
            this.lbExpHigh.Size = new System.Drawing.Size(62, 13);
            this.lbExpHigh.TabIndex = 13;
            this.lbExpHigh.Text = "Exp. High";
            this.lbExpHigh.TextAlign = System.Drawing.ContentAlignment.MiddleCenter;
            // 
            // pnCamControls1
            // 
            this.pnCamControls1.Controls.Add(this.lbNofShots);
            this.pnCamControls1.Controls.Add(this.nmNofShots);
            this.pnCamControls1.Controls.Add(this.sbNofShots);
            this.pnCamControls1.Controls.Add(this.lbExpInc);
            this.pnCamControls1.Controls.Add(this.nmExpInc);
            this.pnCamControls1.Controls.Add(this.sbExpInc);
            this.pnCamControls1.Controls.Add(this.lbExpHigh);
            this.pnCamControls1.Controls.Add(this.sbExpLow);
            this.pnCamControls1.Controls.Add(this.lbExpLow);
            this.pnCamControls1.Controls.Add(this.sbGain);
            this.pnCamControls1.Controls.Add(this.lbExp);
            this.pnCamControls1.Controls.Add(this.nmGain);
            this.pnCamControls1.Controls.Add(this.lbGain);
            this.pnCamControls1.Controls.Add(this.sbExp);
            this.pnCamControls1.Controls.Add(this.nmExpHigh);
            this.pnCamControls1.Controls.Add(this.nmExp);
            this.pnCamControls1.Controls.Add(this.sbExpHigh);
            this.pnCamControls1.Controls.Add(this.nmExpLow);
            this.pnCamControls1.Enabled = false;
            this.pnCamControls1.Location = new System.Drawing.Point(12, 26);
            this.pnCamControls1.Name = "pnCamControls1";
            this.pnCamControls1.Size = new System.Drawing.Size(433, 201);
            this.pnCamControls1.TabIndex = 14;
            // 
            // lbNofShots
            // 
            this.lbNofShots.AutoSize = true;
            this.lbNofShots.Font = new System.Drawing.Font("Microsoft Sans Serif", 8F, System.Drawing.FontStyle.Bold, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.lbNofShots.Location = new System.Drawing.Point(317, 12);
            this.lbNofShots.Name = "lbNofShots";
            this.lbNofShots.Size = new System.Drawing.Size(49, 13);
            this.lbNofShots.TabIndex = 19;
            this.lbNofShots.Text = "# shots";
            this.lbNofShots.TextAlign = System.Drawing.ContentAlignment.MiddleCenter;
            // 
            // nmNofShots
            // 
            this.nmNofShots.Location = new System.Drawing.Point(317, 166);
            this.nmNofShots.Name = "nmNofShots";
            this.nmNofShots.Size = new System.Drawing.Size(50, 20);
            this.nmNofShots.TabIndex = 18;
            this.nmNofShots.ValueChanged += new System.EventHandler(this.nmNofShots_ValueChanged);
            // 
            // sbNofShots
            // 
            this.sbNofShots.Location = new System.Drawing.Point(317, 40);
            this.sbNofShots.Maximum = 35;
            this.sbNofShots.Minimum = 1;
            this.sbNofShots.Name = "sbNofShots";
            this.sbNofShots.Size = new System.Drawing.Size(50, 111);
            this.sbNofShots.TabIndex = 17;
            this.sbNofShots.Value = 16;
            this.sbNofShots.Scroll += new System.Windows.Forms.ScrollEventHandler(this.sbNofShots_Scroll);
            // 
            // lbExpInc
            // 
            this.lbExpInc.AutoSize = true;
            this.lbExpInc.Font = new System.Drawing.Font("Microsoft Sans Serif", 8F, System.Drawing.FontStyle.Bold, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.lbExpInc.Location = new System.Drawing.Point(255, 12);
            this.lbExpInc.Name = "lbExpInc";
            this.lbExpInc.Size = new System.Drawing.Size(54, 13);
            this.lbExpInc.TabIndex = 16;
            this.lbExpInc.Text = "Exp. Inc";
            this.lbExpInc.TextAlign = System.Drawing.ContentAlignment.MiddleCenter;
            this.lbExpInc.Click += new System.EventHandler(this.lbIncGain_Click);
            // 
            // nmExpInc
            // 
            this.nmExpInc.Location = new System.Drawing.Point(257, 166);
            this.nmExpInc.Name = "nmExpInc";
            this.nmExpInc.Size = new System.Drawing.Size(50, 20);
            this.nmExpInc.TabIndex = 15;
            this.nmExpInc.ValueChanged += new System.EventHandler(this.nmExpInc_ValueChanged);
            // 
            // sbExpInc
            // 
            this.sbExpInc.Location = new System.Drawing.Point(257, 40);
            this.sbExpInc.Maximum = 35;
            this.sbExpInc.Minimum = 1;
            this.sbExpInc.Name = "sbExpInc";
            this.sbExpInc.Size = new System.Drawing.Size(50, 111);
            this.sbExpInc.TabIndex = 14;
            this.sbExpInc.Value = 16;
            this.sbExpInc.Scroll += new System.Windows.Forms.ScrollEventHandler(this.sbExpInc_Scroll);
            // 
            // btStart
            // 
            this.btStart.Location = new System.Drawing.Point(704, 465);
            this.btStart.Name = "btStart";
            this.btStart.Size = new System.Drawing.Size(194, 59);
            this.btStart.TabIndex = 15;
            this.btStart.Text = "Start Capturing";
            this.btStart.UseVisualStyleBackColor = true;
            this.btStart.Click += new System.EventHandler(this.btStart_Click);
            // 
            // btStop
            // 
            this.btStop.Location = new System.Drawing.Point(704, 532);
            this.btStop.Name = "btStop";
            this.btStop.Size = new System.Drawing.Size(194, 59);
            this.btStop.TabIndex = 16;
            this.btStop.Text = "Stop Capturing";
            this.btStop.UseVisualStyleBackColor = true;
            this.btStop.Click += new System.EventHandler(this.btStop_Click);
            // 
            // statusStrip1
            // 
            this.statusStrip1.Items.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.tsStatusLabel});
            this.statusStrip1.Location = new System.Drawing.Point(0, 596);
            this.statusStrip1.Name = "statusStrip1";
            this.statusStrip1.Size = new System.Drawing.Size(902, 22);
            this.statusStrip1.TabIndex = 17;
            this.statusStrip1.Text = "statusStrip1";
            // 
            // tsStatusLabel
            // 
            this.tsStatusLabel.Name = "tsStatusLabel";
            this.tsStatusLabel.Size = new System.Drawing.Size(42, 17);
            this.tsStatusLabel.Text = "Status:";
            // 
            // ImageAcquisitionGUI
            // 
            this.AccessibleDescription = "Automated Image Acquisiton Tool for VISMO";
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(902, 618);
            this.Controls.Add(this.statusStrip1);
            this.Controls.Add(this.btStop);
            this.Controls.Add(this.btStart);
            this.Controls.Add(this.pbImage);
            this.Controls.Add(this.pnCamControls1);
            this.Name = "ImageAcquisitionGUI";
            this.Text = "ImageAcquisitionTool";
            this.FormClosing += new System.Windows.Forms.FormClosingEventHandler(this.ImageAcquisitionGUI_FormClosing);
            this.FormClosed += new System.Windows.Forms.FormClosedEventHandler(this.ImageAcquisitionGUI_FormClosed);
            this.Load += new System.EventHandler(this.ImageAcquisitionGUI_Load);
            ((System.ComponentModel.ISupportInitialize)(this.pbImage)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nmGain)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nmExp)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nmExpLow)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nmExpHigh)).EndInit();
            this.pnCamControls1.ResumeLayout(false);
            this.pnCamControls1.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)(this.nmNofShots)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.nmExpInc)).EndInit();
            this.statusStrip1.ResumeLayout(false);
            this.statusStrip1.PerformLayout();
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.PictureBox pbImage;
        private System.Windows.Forms.VScrollBar sbGain;
        private System.Windows.Forms.NumericUpDown nmGain;
        private System.Windows.Forms.NumericUpDown nmExp;
        private System.Windows.Forms.VScrollBar sbExp;
        private System.Windows.Forms.NumericUpDown nmExpLow;
        private System.Windows.Forms.VScrollBar sbExpLow;
        private System.Windows.Forms.NumericUpDown nmExpHigh;
        private System.Windows.Forms.VScrollBar sbExpHigh;
        private System.Windows.Forms.Label lbGain;
        private System.Windows.Forms.Label lbExp;
        private System.Windows.Forms.Label lbExpLow;
        private System.Windows.Forms.Label lbExpHigh;
        private System.Windows.Forms.Panel pnCamControls1;
        private System.Windows.Forms.Label lbNofShots;
        private System.Windows.Forms.NumericUpDown nmNofShots;
        private System.Windows.Forms.VScrollBar sbNofShots;
        private System.Windows.Forms.Label lbExpInc;
        private System.Windows.Forms.NumericUpDown nmExpInc;
        private System.Windows.Forms.VScrollBar sbExpInc;
        private System.Windows.Forms.Button btStart;
        private System.Windows.Forms.Button btStop;
        private System.Windows.Forms.StatusStrip statusStrip1;
        private System.Windows.Forms.ToolStripStatusLabel tsStatusLabel;
    }
}

