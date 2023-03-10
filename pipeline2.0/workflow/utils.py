# -*- coding: utf-8 -*-
# @Author: jsgounot
# @Date:   2023-03-07 13:23:24
# @Last Modified by:   jsgounot
# @Last Modified time: 2023-03-10 09:17:51

import os

def get_db_config(config, name, default):
	try:
		return config['db'][name]
	except:
		return default

def get_compressed_ref(config, sample):
	ref = config['samples'][sample]['ref']
	if ref.endswith('gz'):
		return ref
	else:
		return f'mapping/{sample}/refs/decoys/ref.comp.fa.gz'

def get_ref_path_vcalling(wc, suffix=''):
	refori = 'wodecoys' if wc['bamprefix'].startswith('wodecoys') else 'decoys'
	return f'mapping/{wc.sample}/refs/{refori}/ref.fa.bgz' + suffix

def get_target_r1(config, sample):
	target = config['samples'][sample]['target']
	return f'simulation/illumina/{sample}/{target}_R1.fq'

def get_target_r2(config, sample):
	target = config['samples'][sample]['target']
	return f'simulation/illumina/{sample}/{target}_R2.fq'

def get_target_mutated_ref(wc, config):
	target = config['samples'][wc.sample]['target']
	return f'simulation/mutated/{wc.sample}/{target}.fa'