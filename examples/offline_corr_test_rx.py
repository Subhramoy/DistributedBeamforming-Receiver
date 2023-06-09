#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Offline Corr Test Rx
# Generated: Thu Oct 31 23:25:39 2019
##################################################

if __name__ == '__main__':
    import ctypes
    import sys
    if sys.platform.startswith('linux'):
        try:
            x11 = ctypes.cdll.LoadLibrary('libX11.so')
            x11.XInitThreads()
        except:
            print "Warning: failed to XInitThreads()"

from PyQt4 import Qt
from gnuradio import analog
from gnuradio import blocks
from gnuradio import eng_notation
from gnuradio import gr
from gnuradio import qtgui
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from optparse import OptionParser
import beamforming
import sip
import sys
from gnuradio import qtgui


class offline_corr_test_rx(gr.top_block, Qt.QWidget):

    def __init__(self):
        gr.top_block.__init__(self, "Offline Corr Test Rx")
        Qt.QWidget.__init__(self)
        self.setWindowTitle("Offline Corr Test Rx")
        qtgui.util.check_set_qss()
        try:
            self.setWindowIcon(Qt.QIcon.fromTheme('gnuradio-grc'))
        except:
            pass
        self.top_scroll_layout = Qt.QVBoxLayout()
        self.setLayout(self.top_scroll_layout)
        self.top_scroll = Qt.QScrollArea()
        self.top_scroll.setFrameStyle(Qt.QFrame.NoFrame)
        self.top_scroll_layout.addWidget(self.top_scroll)
        self.top_scroll.setWidgetResizable(True)
        self.top_widget = Qt.QWidget()
        self.top_scroll.setWidget(self.top_widget)
        self.top_layout = Qt.QVBoxLayout(self.top_widget)
        self.top_grid_layout = Qt.QGridLayout()
        self.top_layout.addLayout(self.top_grid_layout)

        self.settings = Qt.QSettings("GNU Radio", "offline_corr_test_rx")
        self.restoreGeometry(self.settings.value("geometry").toByteArray())


        ##################################################
        # Variables
        ##################################################
        self.trainingSignal_size = trainingSignal_size = 16456
        self.samp_rate = samp_rate = 32000
        self.data_files_path = data_files_path = "/home/genesys/workarea-gnuradio/gnuradio/gr-beamforming/examples/data"

        ##################################################
        # Blocks
        ##################################################
        self.zero_padding_0_0_0_1 = analog.sig_source_c(0, analog.GR_CONST_WAVE, 0, 0, 0.1-0.1j)
        self.zero_padding_0_0_0_0_1 = analog.sig_source_c(0, analog.GR_CONST_WAVE, 0, 0, 0)
        self.zero_padding_0_0_0_0_0 = analog.sig_source_c(0, analog.GR_CONST_WAVE, 0, 0, 0)
        self.zero_padding_0_0_0_0 = analog.sig_source_c(0, analog.GR_CONST_WAVE, 0, 0, 0)
        self.zero_padding_0_0_0 = analog.sig_source_c(0, analog.GR_CONST_WAVE, 0, 0, 0.1-0.1j)
        self.qtgui_time_sink_x_0_0_1_1_0_0_1 = qtgui.time_sink_c(
        	200000, #size
        	samp_rate, #samp_rate
        	'Retrieved Payload', #name
        	1 #number of inputs
        )
        self.qtgui_time_sink_x_0_0_1_1_0_0_1.set_update_time(0.10)
        self.qtgui_time_sink_x_0_0_1_1_0_0_1.set_y_axis(-1, 1)

        self.qtgui_time_sink_x_0_0_1_1_0_0_1.set_y_label('Amplitude', "")

        self.qtgui_time_sink_x_0_0_1_1_0_0_1.enable_tags(-1, True)
        self.qtgui_time_sink_x_0_0_1_1_0_0_1.set_trigger_mode(qtgui.TRIG_MODE_FREE, qtgui.TRIG_SLOPE_POS, 0.0, 0, 0, "")
        self.qtgui_time_sink_x_0_0_1_1_0_0_1.enable_autoscale(True)
        self.qtgui_time_sink_x_0_0_1_1_0_0_1.enable_grid(False)
        self.qtgui_time_sink_x_0_0_1_1_0_0_1.enable_axis_labels(True)
        self.qtgui_time_sink_x_0_0_1_1_0_0_1.enable_control_panel(False)
        self.qtgui_time_sink_x_0_0_1_1_0_0_1.enable_stem_plot(False)

        if not True:
          self.qtgui_time_sink_x_0_0_1_1_0_0_1.disable_legend()

        labels = ['IQ', 'Corr Output', '', '', '',
                  '', '', '', '', '']
        widths = [1, 1, 1, 1, 1,
                  1, 1, 1, 1, 1]
        colors = ["blue", "red", "green", "black", "cyan",
                  "magenta", "yellow", "dark red", "dark green", "blue"]
        styles = [1, 1, 1, 1, 1,
                  1, 1, 1, 1, 1]
        markers = [-1, -1, -1, -1, -1,
                   -1, -1, -1, -1, -1]
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
                  1.0, 1.0, 1.0, 1.0, 1.0]

        for i in xrange(2):
            if len(labels[i]) == 0:
                if(i % 2 == 0):
                    self.qtgui_time_sink_x_0_0_1_1_0_0_1.set_line_label(i, "Re{{Data {0}}}".format(i/2))
                else:
                    self.qtgui_time_sink_x_0_0_1_1_0_0_1.set_line_label(i, "Im{{Data {0}}}".format(i/2))
            else:
                self.qtgui_time_sink_x_0_0_1_1_0_0_1.set_line_label(i, labels[i])
            self.qtgui_time_sink_x_0_0_1_1_0_0_1.set_line_width(i, widths[i])
            self.qtgui_time_sink_x_0_0_1_1_0_0_1.set_line_color(i, colors[i])
            self.qtgui_time_sink_x_0_0_1_1_0_0_1.set_line_style(i, styles[i])
            self.qtgui_time_sink_x_0_0_1_1_0_0_1.set_line_marker(i, markers[i])
            self.qtgui_time_sink_x_0_0_1_1_0_0_1.set_line_alpha(i, alphas[i])

        self._qtgui_time_sink_x_0_0_1_1_0_0_1_win = sip.wrapinstance(self.qtgui_time_sink_x_0_0_1_1_0_0_1.pyqwidget(), Qt.QWidget)
        self.top_grid_layout.addWidget(self._qtgui_time_sink_x_0_0_1_1_0_0_1_win)
        self.qtgui_time_sink_x_0_0_1_1_0_0_0 = qtgui.time_sink_f(
        	200000, #size
        	samp_rate, #samp_rate
        	'XCor', #name
        	1 #number of inputs
        )
        self.qtgui_time_sink_x_0_0_1_1_0_0_0.set_update_time(0.10)
        self.qtgui_time_sink_x_0_0_1_1_0_0_0.set_y_axis(-1, 1)

        self.qtgui_time_sink_x_0_0_1_1_0_0_0.set_y_label('Amplitude', "")

        self.qtgui_time_sink_x_0_0_1_1_0_0_0.enable_tags(-1, True)
        self.qtgui_time_sink_x_0_0_1_1_0_0_0.set_trigger_mode(qtgui.TRIG_MODE_FREE, qtgui.TRIG_SLOPE_POS, 0.0, 0, 0, "")
        self.qtgui_time_sink_x_0_0_1_1_0_0_0.enable_autoscale(True)
        self.qtgui_time_sink_x_0_0_1_1_0_0_0.enable_grid(False)
        self.qtgui_time_sink_x_0_0_1_1_0_0_0.enable_axis_labels(True)
        self.qtgui_time_sink_x_0_0_1_1_0_0_0.enable_control_panel(False)
        self.qtgui_time_sink_x_0_0_1_1_0_0_0.enable_stem_plot(False)

        if not True:
          self.qtgui_time_sink_x_0_0_1_1_0_0_0.disable_legend()

        labels = ['IQ', 'Corr Output', '', '', '',
                  '', '', '', '', '']
        widths = [1, 1, 1, 1, 1,
                  1, 1, 1, 1, 1]
        colors = ["blue", "red", "green", "black", "cyan",
                  "magenta", "yellow", "dark red", "dark green", "blue"]
        styles = [1, 1, 1, 1, 1,
                  1, 1, 1, 1, 1]
        markers = [-1, -1, -1, -1, -1,
                   -1, -1, -1, -1, -1]
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
                  1.0, 1.0, 1.0, 1.0, 1.0]

        for i in xrange(1):
            if len(labels[i]) == 0:
                self.qtgui_time_sink_x_0_0_1_1_0_0_0.set_line_label(i, "Data {0}".format(i))
            else:
                self.qtgui_time_sink_x_0_0_1_1_0_0_0.set_line_label(i, labels[i])
            self.qtgui_time_sink_x_0_0_1_1_0_0_0.set_line_width(i, widths[i])
            self.qtgui_time_sink_x_0_0_1_1_0_0_0.set_line_color(i, colors[i])
            self.qtgui_time_sink_x_0_0_1_1_0_0_0.set_line_style(i, styles[i])
            self.qtgui_time_sink_x_0_0_1_1_0_0_0.set_line_marker(i, markers[i])
            self.qtgui_time_sink_x_0_0_1_1_0_0_0.set_line_alpha(i, alphas[i])

        self._qtgui_time_sink_x_0_0_1_1_0_0_0_win = sip.wrapinstance(self.qtgui_time_sink_x_0_0_1_1_0_0_0.pyqwidget(), Qt.QWidget)
        self.top_grid_layout.addWidget(self._qtgui_time_sink_x_0_0_1_1_0_0_0_win)
        self.qtgui_time_sink_x_0_0_1_1_0_0 = qtgui.time_sink_c(
        	200000, #size
        	samp_rate, #samp_rate
        	'Rx Signal', #name
        	1 #number of inputs
        )
        self.qtgui_time_sink_x_0_0_1_1_0_0.set_update_time(0.10)
        self.qtgui_time_sink_x_0_0_1_1_0_0.set_y_axis(-1, 1)

        self.qtgui_time_sink_x_0_0_1_1_0_0.set_y_label('Amplitude', "")

        self.qtgui_time_sink_x_0_0_1_1_0_0.enable_tags(-1, True)
        self.qtgui_time_sink_x_0_0_1_1_0_0.set_trigger_mode(qtgui.TRIG_MODE_FREE, qtgui.TRIG_SLOPE_POS, 0.0, 0, 0, "")
        self.qtgui_time_sink_x_0_0_1_1_0_0.enable_autoscale(True)
        self.qtgui_time_sink_x_0_0_1_1_0_0.enable_grid(False)
        self.qtgui_time_sink_x_0_0_1_1_0_0.enable_axis_labels(True)
        self.qtgui_time_sink_x_0_0_1_1_0_0.enable_control_panel(False)
        self.qtgui_time_sink_x_0_0_1_1_0_0.enable_stem_plot(False)

        if not True:
          self.qtgui_time_sink_x_0_0_1_1_0_0.disable_legend()

        labels = ['IQ', 'Corr Output', '', '', '',
                  '', '', '', '', '']
        widths = [1, 1, 1, 1, 1,
                  1, 1, 1, 1, 1]
        colors = ["blue", "red", "green", "black", "cyan",
                  "magenta", "yellow", "dark red", "dark green", "blue"]
        styles = [1, 1, 1, 1, 1,
                  1, 1, 1, 1, 1]
        markers = [-1, -1, -1, -1, -1,
                   -1, -1, -1, -1, -1]
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
                  1.0, 1.0, 1.0, 1.0, 1.0]

        for i in xrange(2):
            if len(labels[i]) == 0:
                if(i % 2 == 0):
                    self.qtgui_time_sink_x_0_0_1_1_0_0.set_line_label(i, "Re{{Data {0}}}".format(i/2))
                else:
                    self.qtgui_time_sink_x_0_0_1_1_0_0.set_line_label(i, "Im{{Data {0}}}".format(i/2))
            else:
                self.qtgui_time_sink_x_0_0_1_1_0_0.set_line_label(i, labels[i])
            self.qtgui_time_sink_x_0_0_1_1_0_0.set_line_width(i, widths[i])
            self.qtgui_time_sink_x_0_0_1_1_0_0.set_line_color(i, colors[i])
            self.qtgui_time_sink_x_0_0_1_1_0_0.set_line_style(i, styles[i])
            self.qtgui_time_sink_x_0_0_1_1_0_0.set_line_marker(i, markers[i])
            self.qtgui_time_sink_x_0_0_1_1_0_0.set_line_alpha(i, alphas[i])

        self._qtgui_time_sink_x_0_0_1_1_0_0_win = sip.wrapinstance(self.qtgui_time_sink_x_0_0_1_1_0_0.pyqwidget(), Qt.QWidget)
        self.top_grid_layout.addWidget(self._qtgui_time_sink_x_0_0_1_1_0_0_win)
        self.qtgui_freq_sink_x_0_0_0_0 = qtgui.freq_sink_c(
        	256*8, #size
        	firdes.WIN_BLACKMAN_hARRIS, #wintype
        	900e6, #fc
        	samp_rate, #bw
        	'Rx Signal', #name
        	1 #number of inputs
        )
        self.qtgui_freq_sink_x_0_0_0_0.set_update_time(0.10)
        self.qtgui_freq_sink_x_0_0_0_0.set_y_axis(-140, 10)
        self.qtgui_freq_sink_x_0_0_0_0.set_y_label('Relative Gain', 'dB')
        self.qtgui_freq_sink_x_0_0_0_0.set_trigger_mode(qtgui.TRIG_MODE_FREE, 0.0, 0, "")
        self.qtgui_freq_sink_x_0_0_0_0.enable_autoscale(False)
        self.qtgui_freq_sink_x_0_0_0_0.enable_grid(False)
        self.qtgui_freq_sink_x_0_0_0_0.set_fft_average(1.0)
        self.qtgui_freq_sink_x_0_0_0_0.enable_axis_labels(True)
        self.qtgui_freq_sink_x_0_0_0_0.enable_control_panel(False)

        if not True:
          self.qtgui_freq_sink_x_0_0_0_0.disable_legend()

        if "complex" == "float" or "complex" == "msg_float":
          self.qtgui_freq_sink_x_0_0_0_0.set_plot_pos_half(not True)

        labels = ['', '', '', '', '',
                  '', '', '', '', '']
        widths = [1, 1, 1, 1, 1,
                  1, 1, 1, 1, 1]
        colors = ["blue", "red", "green", "black", "cyan",
                  "magenta", "yellow", "dark red", "dark green", "dark blue"]
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
                  1.0, 1.0, 1.0, 1.0, 1.0]
        for i in xrange(1):
            if len(labels[i]) == 0:
                self.qtgui_freq_sink_x_0_0_0_0.set_line_label(i, "Data {0}".format(i))
            else:
                self.qtgui_freq_sink_x_0_0_0_0.set_line_label(i, labels[i])
            self.qtgui_freq_sink_x_0_0_0_0.set_line_width(i, widths[i])
            self.qtgui_freq_sink_x_0_0_0_0.set_line_color(i, colors[i])
            self.qtgui_freq_sink_x_0_0_0_0.set_line_alpha(i, alphas[i])

        self._qtgui_freq_sink_x_0_0_0_0_win = sip.wrapinstance(self.qtgui_freq_sink_x_0_0_0_0.pyqwidget(), Qt.QWidget)
        self.top_grid_layout.addWidget(self._qtgui_freq_sink_x_0_0_0_0_win)
        self.blocks_throttle_0 = blocks.throttle(gr.sizeof_gr_complex*1, samp_rate,True)
        self.blocks_tag_gate_0 = blocks.tag_gate(gr.sizeof_gr_complex * 1, False)
        self.blocks_tag_gate_0.set_single_key("")
        self.blocks_stream_mux_1_1 = blocks.stream_mux(gr.sizeof_gr_complex*1, (trainingSignal_size, 400 , 256* 64 , 100))
        self.blocks_stream_mux_1_0 = blocks.stream_mux(gr.sizeof_gr_complex*1, (1000 , 256* 64))
        self.blocks_stream_mux_1 = blocks.stream_mux(gr.sizeof_gr_complex*1, (trainingSignal_size, 400 , 256* 64 , 100))
        self.blocks_repeat_0_2 = blocks.repeat(gr.sizeof_gr_complex*1, 400)
        self.blocks_repeat_0_1_0 = blocks.repeat(gr.sizeof_gr_complex*1, 100)
        self.blocks_repeat_0_1 = blocks.repeat(gr.sizeof_gr_complex*1, 100)
        self.blocks_repeat_0_0_0 = blocks.repeat(gr.sizeof_gr_complex*1, 256*64)
        self.blocks_repeat_0_0 = blocks.repeat(gr.sizeof_gr_complex*1, 256*64)
        self.blocks_repeat_0 = blocks.repeat(gr.sizeof_gr_complex*1, 400)
        self.blocks_complex_to_mag_0 = blocks.complex_to_mag(1)
        self.blocks_add_xx_0 = blocks.add_vcc(1)
        self.beamforming_matlab_file_payload_py_0_0 = beamforming.matlab_file_payload_py(data_files_path + "/trainingSig2")
        self.beamforming_matlab_file_payload_py_0 = beamforming.matlab_file_payload_py(data_files_path + "/trainingSig1")
        self.beamforming_filter_payload_py_0 = beamforming.filter_payload_py('payload')
        self.beamforming_correlate_and_tag_py_0 = beamforming.correlate_and_tag_py(trainingSignal_size, trainingSignal_size + 400 + 256* 64 + 100, 2, data_files_path + "/trainingSig", 1, 0)



        ##################################################
        # Connections
        ##################################################
        self.connect((self.beamforming_correlate_and_tag_py_0, 0), (self.beamforming_filter_payload_py_0, 0))
        self.connect((self.beamforming_correlate_and_tag_py_0, 1), (self.blocks_complex_to_mag_0, 0))
        self.connect((self.beamforming_correlate_and_tag_py_0, 0), (self.qtgui_freq_sink_x_0_0_0_0, 0))
        self.connect((self.beamforming_correlate_and_tag_py_0, 0), (self.qtgui_time_sink_x_0_0_1_1_0_0, 0))
        self.connect((self.beamforming_filter_payload_py_0, 0), (self.blocks_stream_mux_1_0, 1))
        self.connect((self.beamforming_matlab_file_payload_py_0, 0), (self.blocks_stream_mux_1, 0))
        self.connect((self.beamforming_matlab_file_payload_py_0_0, 0), (self.blocks_stream_mux_1_1, 0))
        self.connect((self.blocks_add_xx_0, 0), (self.blocks_throttle_0, 0))
        self.connect((self.blocks_complex_to_mag_0, 0), (self.qtgui_time_sink_x_0_0_1_1_0_0_0, 0))
        self.connect((self.blocks_repeat_0, 0), (self.blocks_stream_mux_1, 1))
        self.connect((self.blocks_repeat_0_0, 0), (self.blocks_stream_mux_1, 2))
        self.connect((self.blocks_repeat_0_0_0, 0), (self.blocks_stream_mux_1_1, 2))
        self.connect((self.blocks_repeat_0_1, 0), (self.blocks_stream_mux_1, 3))
        self.connect((self.blocks_repeat_0_1_0, 0), (self.blocks_stream_mux_1_1, 3))
        self.connect((self.blocks_repeat_0_2, 0), (self.blocks_stream_mux_1_1, 1))
        self.connect((self.blocks_stream_mux_1, 0), (self.blocks_add_xx_0, 0))
        self.connect((self.blocks_stream_mux_1_0, 0), (self.blocks_tag_gate_0, 0))
        self.connect((self.blocks_stream_mux_1_1, 0), (self.blocks_add_xx_0, 1))
        self.connect((self.blocks_tag_gate_0, 0), (self.qtgui_time_sink_x_0_0_1_1_0_0_1, 0))
        self.connect((self.blocks_throttle_0, 0), (self.beamforming_correlate_and_tag_py_0, 0))
        self.connect((self.zero_padding_0_0_0, 0), (self.blocks_repeat_0_0, 0))
        self.connect((self.zero_padding_0_0_0_0, 0), (self.blocks_repeat_0, 0))
        self.connect((self.zero_padding_0_0_0_0, 0), (self.blocks_repeat_0_1, 0))
        self.connect((self.zero_padding_0_0_0_0_0, 0), (self.blocks_stream_mux_1_0, 0))
        self.connect((self.zero_padding_0_0_0_0_1, 0), (self.blocks_repeat_0_1_0, 0))
        self.connect((self.zero_padding_0_0_0_0_1, 0), (self.blocks_repeat_0_2, 0))
        self.connect((self.zero_padding_0_0_0_1, 0), (self.blocks_repeat_0_0_0, 0))

    def closeEvent(self, event):
        self.settings = Qt.QSettings("GNU Radio", "offline_corr_test_rx")
        self.settings.setValue("geometry", self.saveGeometry())
        event.accept()

    def get_trainingSignal_size(self):
        return self.trainingSignal_size

    def set_trainingSignal_size(self, trainingSignal_size):
        self.trainingSignal_size = trainingSignal_size

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.qtgui_time_sink_x_0_0_1_1_0_0_1.set_samp_rate(self.samp_rate)
        self.qtgui_time_sink_x_0_0_1_1_0_0_0.set_samp_rate(self.samp_rate)
        self.qtgui_time_sink_x_0_0_1_1_0_0.set_samp_rate(self.samp_rate)
        self.qtgui_freq_sink_x_0_0_0_0.set_frequency_range(900e6, self.samp_rate)
        self.blocks_throttle_0.set_sample_rate(self.samp_rate)

    def get_data_files_path(self):
        return self.data_files_path

    def set_data_files_path(self, data_files_path):
        self.data_files_path = data_files_path


def main(top_block_cls=offline_corr_test_rx, options=None):

    from distutils.version import StrictVersion
    if StrictVersion(Qt.qVersion()) >= StrictVersion("4.5.0"):
        style = gr.prefs().get_string('qtgui', 'style', 'raster')
        Qt.QApplication.setGraphicsSystem(style)
    qapp = Qt.QApplication(sys.argv)

    tb = top_block_cls()
    tb.start()
    tb.show()

    def quitting():
        tb.stop()
        tb.wait()
    qapp.connect(qapp, Qt.SIGNAL("aboutToQuit()"), quitting)
    qapp.exec_()


if __name__ == '__main__':
    main()
