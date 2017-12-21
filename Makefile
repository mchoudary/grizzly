# Makefile for code and data in grizzly project (template attacks on 8-bit XMEGA 256 A3U)
# See the readme.txt file for more details.
#
# Author: Omar Choudary (omar.choudary@cl.cam.ac.uk)

FDATA=e2_bat_fb_beta_raw_s_0_3071.raw
MATLAB=matlab -nosplash

all: figures/a2_bat_fb_dlinear_n200r_ls_r10_guess_entropy.pdf

data: $(FDATA)

$(FDATA):
	[ -e $@ ] || \
	wget http://www.cl.cam.ac.uk/research/security/datasets/grizzly/$(FDATA).gz; \
	gunzip $(FDATA).gz

results/a2_bat_fb_templates_dlinear_n200r_slr_g1000_r10.mat: \
		do_test_success_templates_bat_fb_dlinear.m \
		$(FDATA)
	[ -d results ] || mkdir results
	echo do_test_success_templates_bat_fb_dlinear | $(MATLAB)

figures/a2_bat_fb_dlinear_n200r_ls_r10_guess_entropy.pdf: \
		do_show_results_templates_a2_bat_fb.m \
		results/a2_bat_fb_templates_dlinear_n200r_slr_g1000_r10.mat
	[ -d figures ] || mkdir figures
	echo do_show_results_templates_a2_bat_fb | $(MATLAB)

clean:
	rm -f figures/a2_bat_fb_dlinear_n200r_ls_r10_guess_entropy.pdf
