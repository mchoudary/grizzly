DATA and CODE for project GRIZZLY: efficient and portable template attacks.
Author: Marios Omar Choudary (omar.choudary@cl.cam.ac.uk,marios.choudary@cs.pub.ro)
Version: 1.2
Last updated: 22 December 2017

PAPERS

The data and code in this folder allow you to reproduce the results shown
in these papers:

1) Marios O. Choudary and Markus G. Kuhn, "Efficient, Portable Template Attacks",
IEEE Transactions on Information Forensics and Security, 13(2), 2018, pp.490--501.

2) Omar Choudary and Markus G. Kuhn, "Template attacks on different devices",
COSADE 2014, Paris, 2014. Springer LNCS 8622, pp. 179--198.

3) Omar Choudary and Markus G. Kuhn, "Efficient Template Attacks", CARDIS 2013,
Berlin, 27--29 November, 2013. Springer LNCS 8419, pp. 253--270.

If you find the associated data or code useful please cite the relevant papers
in your results.

DATA

You can download the data associated with these papers here:
https://www.cl.cam.ac.uk/research/security/datasets/grizzly/

Then you should uncompress each of the data files by running gunzip, e.g.:
"gunzip e2_bat_fb_beta_raw_s_0_3071.raw.gz"

and put the resulting file (e.g. e2_bat_fb_beta_raw_s_0_3071.raw) in the same folder
as this readme file (i.e. the folder with all the MATLAB scripts).

Each of the uncompressed data files contain the raw traces that we
acquired from our four evaluation boards (Alpha, Beta, Gamma, Delta) on the
Atmel XMEGA A3U microcontroller, by using a Tektronix TDS 7054 oscilloscope
using the full bandwidth (500MHz) and batteries via a voltage regulator 3.3 V
as power supply. There are five datasets available: four corresponding to
each device (alpha, beta, gamma, delta) and a second dataset on device Beta,
which corresponds to a different acquisition campaign to observe similarities
and differences between attacks on different devices and attacks on different
acquisition campaigns on the same device (see papers for details).

The full acquisition settings for all campaigns/devices are:
Power supply from batteries via voltage regulator at 3.3V,
1MHz sine clk (3 Vpp,1.5 V DC offset, 50 Ohm load, sine),
10 Ohm resistor, active probe. Bandwidth 500MHz,250MS/s,4ns/point,1us/div,
2500 pt/frame,SAMPLE mode with Fastframe, 10mV/div vertical resolution,
DC coupling.

The experiments (which we called E2) are as follows:
Loading bytes into registers r8:r15. All bytes are 0 except for the
second byte that gets loaded into r9, which is changed between 0 and 255.
The experiment is as follows:
  - for nr_trials=1310720=256*5120 times do:
    1. initialize registers with zeros.
    2. copy 8 bytes of data into r8:r15 at address A=0x2525
  - each block of 8 bytes has first byte 0, next byte is a value between 0
    and 255 and the next 6 bytes are 0.
  - every 256 trials all the 256 values {0, ..., 255} are used for byte 2.
    That is, a random permutation of the values {0, ..., 255} is used as
    byte 2 over 256 trials, so we use every byte value for 5120 times.

This is the code that is being run from the point of trigger:
// signal start of loading data trigger (PC1)
sbi _SFR_IO_ADDR(VPORT0_DIR), _SFR_IO_ADDR(PIN1_bp)     ; 1 clock cycle
sbi _SFR_IO_ADDR(VPORT0_OUT), _SFR_IO_ADDR(PIN1_bp)     ; 1 clock cycle
nop
cbi _SFR_IO_ADDR(VPORT0_OUT), _SFR_IO_ADDR(PIN1_bp)     ; 1 clock cycle
cbi _SFR_IO_ADDR(VPORT0_DIR), _SFR_IO_ADDR(PIN1_bp)     ; 1 clock cycle
// add a few more nops after signal to remove influence on power consumption
// there are 6 clock cycles between trigger and mov instruction
nop
nop
nop

// load the data
movw r30, r24   ; 1 clock cycle
ld r8, Z+       ; 2 clock cycles per ld when accessing internal RAM
ld r9, Z+
ld r10, Z+
ld r11, Z+
ld r12, Z+
ld r13, Z+
ld r14, Z+
ld r15, Z+

// 3 nops
nop
nop
nop

The initial values of the registers r8-r15 are 0x00.

In each trace, the point at 5% corresponds to the MOV instruction before the
sequence of LD instructions that load the DES key. Therefore, the point at 15%
corresponds to the first clock cycle of the first LD instruction, the point at
25% to the second clock cycle, the point at 35% to the first clock of the second
LD instruction and so on.

MATLAB SCRIPTS

The MATLAB scripts in this folder will allow you to run template attacks
on the raw data using different parameters. The scripts allow you to run
attacks both on a single device (i.e. data from a single board) as well
as attacks that use profiling from one device but target a different device.

Start by looking at:
'do_test_success_templates_bat_fb_dlinear.m'

which will run the template attacks for a single device (Beta).
The results of this method will be stored in:
'results/a2_bat_fb_templates_dlinear_n200r_slr_g1000_r10.mat'

(you must create the folder "results/". See the Makefile below)

You can then use the script:
'do_show_results_templates_a2_bat_fb.m'

on the results above to produce a figure similar to the ones in our papers.
The figure will be stored in:
'figures/a2_bat_fb_dlinear_n200r_ls_r10_guess_entropy.pdf'

(you must create the folder "figures/". See the Makefile below)

For results using different devices, you should look at the scripts:
'do_test_success_templates_a2d_ab_bat_fb_dlinear.m'
(normal template attacks)

and:
'do_test_success_templates_a2d_ab_adapt_boffset.m'
(template attack using our method to insert some random offset during training)

The associated scripts to produce figures from those attacks are:
'do_show_results_templates_a2d_ab_bat_fb.m'

and:
'do_show_results_templates_a2d_ab_adapt_boffset.m'

These scripts should allow you to reproduce some of the figures in our papers
and they should also allow you to easily build your own scripts to reproduce
all figures (e.g. using also the other campaigns, not just alpha and beta)
as well as to test your own algorithms on our datasets.

AUTOMATIC BUILD

In this folder there should also be a Makefile that allows you to
automatically build things (you need make, MATLAB in your path and perhaps
a Linux environment since I have not tested in any other environment).

Run one of the following:

"make" or "make all":
  download two data sets (alpha and beta), run the template attacks on
  single device (beta) as well as multiple devices (alpha, beta) and
  produce result figures.

"make multi":
  run only the attacks on different devices.

"make single":
  run only the attacks on a single device.

"make data":
  download and uncompress the data only

See the file Makefile for more details.

Enjoy!
