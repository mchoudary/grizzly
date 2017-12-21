DATA and CODE for project GRIZZLY: efficient template attacks.
Author: Omar Choudary (omar.choudary@cl.cam.ac.uk)
Version: 1.1
Last updated: 2 January 2014

PAPER

The data and code in this folder allow you to reproduce the results shown
in this paper:
Omar Choudary and Markus G. Kuhn, "Efficient Template Attacks", CARDIS 2013,
Berlin, 27--29 November, 2013.

If you find the associated data or code useful please cite the above paper
in your results.

DATA

You can download the data associated with this paper from here:
http://www.cl.cam.ac.uk/research/security/datasets/grizzly/e2_bat_fb_beta_raw_s_0_3071.raw.gz

Then you should uncompress this file by running gunzip:
"gunzip e2_bat_fb_beta_raw_s_0_3071.raw.gz"

and put the resulting file (e2_bat_fb_beta_raw_s_0_3071.raw) in the same folder
as this readme file (i.e. the folder with all the MATLAB scripts).

The data in e2_bat_fb_beta_raw_s_0_3071.raw contains the raw traces that we
acquired from our evaluation board on the Atmel XMEGA A3U microcontroller,
by using a Tektronix TDS 7054 oscilloscope using the full bandwidth
(500MHz) and batteries via a voltage regulator 3.3 V as power supply.

The full settings are:
Power supply from batteries via voltage regulator at 3.3V,
1MHz sine clk (3 Vpp,1.5 V DC offset, 50 Ohm load, sine),
10 Ohm resistor, active probe. Bandwidth 500MHz,250MS/s,4ns/point,1us/div,
2500 pt/frame,SAMPLE mode with Fastframe, 10mV/div vertical resolution,
DC coupling.

The experiments (E2) are as follows:
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
on the raw data using different parameters.

Start by looking at do_test_success_templates_bat_fb_dlinear.m, which will
run the template attacks with all the parameters described in the above
paper. The results of this method will be stored in
results/a2_bat_fb_templates_dlinear_n200r_slr_g1000_r10.mat
(you must create the folder "results/". See the Makefile below)

You can then use the script do_show_results_templates_a2_bat_fb.m
on these results to produce a figure similar to the one in our paper.
The figure will be stored in
figures/a2_bat_fb_dlinear_n200r_ls_r10_guess_entropy.pdf
(you must create the folder "figures/". See the Makefile below)

AUTOMATIC BUILD

In this folder there should also be a Makefile that allows you to
automatically build things (you need make, MATLAB in your path and perhaps
a linux environment since I have not tested in any other environment).

Run one of the following:

"make" or "make all": downloads the data, runs the template attacks and produces results figure.
"make data": downloads and uncompresses the data only

See the file Makefile for more details.

Enjoy!
