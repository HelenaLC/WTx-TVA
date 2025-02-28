# WTx-CosMx SMI on TVA/CRC

import re
import pandas as pd

ex = "meta/exclude.rds" # barcodes overlapping with CARD8 barcode in 1st run
md = "data/raw/metadata.txt" # slide metadata (including FOV-to-section mapping)
MD = pd.read_csv(md, dtype="string", delimiter="\t")

DID = []
for i in range(len(MD)):
	did = MD["did"][i]
	if did not in DID:
		DID += [did]
SID = sorted(set(MD["sid"]))
SID = [x for x in SID if x != "241"]
SUB = ["epi", "imm", "str"]

# targets
raw = "outs/raw-{sid}"
fil = "outs/fil-{sid}.rds"
pro = "outs/pro-{sid}.rds"
pol = "outs/pol-{sid}.parquet"
roi = "outs/roi-{sid}.rds"
ccc = "outs/ccc-{sid}.rds"
sig = "outs/sig-{sid}.rds"
#pro = "outs/pro-{sid}.rds"

ROI = {}
for sid in SID:
	did = MD["did"][MD["sid"] == sid]
	foo = "imgs/{did}/shapes/{sid}/{{roi}}.pickle"
	foo = expand(foo, did=did, sid=sid)[0]
	foo = glob_wildcards(foo).roi
	if sid not in ROI.keys():
		ROI[sid] = {}
	ROI[sid] = foo

pbs_lv1 = "data/ref/outs/pbs-lv1.rds"
pbs_lv2 = "data/ref/outs/pbs-lv2,{sub}.rds"

ist = "outs/ist-{sid}.rds"
lv1 = "outs/lv1-{sid}.rds"

sub = "outs/sub-{sid},{sub}.rds"
rep = "outs/rep-{sid},{sub}.rds"
trj = "outs/trj-{sid},{sub}.rds"
add = "outs/add-{sub}"

jst = "outs/jst-{sid},{sub}.rds"
lv2 = "outs/lv2-{sid},{sub}.rds"
pbs = "outs/pbs-{sid},{sub}.rds"
ctx = "outs/ctx-{sid}.rds"

add_by_sid = "outs/add_by_sid-{sid}.rds"
add_by_sub = "outs/add_by_sub-{sub}.rds"

# plotting
plt_raw_all = []
plt_any_all = []
plt_raw = "plts/raw/raw,{out2}-{plt},{sid}.pdf"
plt_raw_ = "plts/raw/raw,{out2}-{plt},{{sid}}.pdf"
plt_any = "plts/any/{out1},{out2}-{plt},{sid}.pdf"
plt_any_ = "plts/any/{out1},{out2}-{plt},{{sid}}.pdf"
foo = glob_wildcards("code/10-plt_{x},{y}-{z}.R")
for x,y,z in zip(foo.x, foo.y, foo.z):
	if x == "raw":
		ps = expand(plt_raw_, out1=x, out2=y, plt=z)
		plt_raw_all += expand(ps, sid=SID)
	elif x != "all":
		ps = expand(plt_any_, out1=x, out2=y, plt=z)
		plt_any_all += expand(ps, sid=SID)

plt = []

plt__sid = "plts/{out1},{plt},{sid}.pdf"
plt__sid_ = "plts/{out1},{plt},{{sid}}.pdf"
foo = glob_wildcards("code/10-plt__sid-{x},{y}.R")
for x,y in zip(foo.x, foo.y):
	plt += expand(plt__sid_, out1=x, plt=y)

plt__sid__sid = "plts/{out1},{out2},{plt},{sid}.pdf"
plt__sid__sid_ = "plts/{out1},{out2},{plt},{{sid}}.pdf"
foo = glob_wildcards("code/10-plt__sid__sid-{x},{y},{z}.R")
for x,y,z in zip(foo.x, foo.y, foo.z):
	plt += expand(plt__sid__sid_, out1=x, out2=y, plt=z)

plt__all_sid = "plts/{out1},{plt}.pdf"
foo = glob_wildcards("code/10-plt__all_sid-{x},{y}.R")
for x,y in zip(foo.x, foo.y):
	plt += expand(plt__all_sid, out1=x, plt=y)

plt__all_sid__all_sid = "plts/{out1},{out2},{plt}.pdf"
foo = glob_wildcards("code/10-plt__all_sid__all_sid-{x},{y},{z}.R")
for x,y,z in zip(foo.x, foo.y, foo.z):
	plt += expand(plt__all_sid__all_sid, out1=x, out2=y, plt=z)

plt__sid_sub = "plts/{out1},{plt},{sid},{sub}.pdf"
plt__sid_sub_ = "plts/{out1},{plt},{{sid}},{{sub}}.pdf"
foo = glob_wildcards("code/10-plt__sid_sub-{x},{y}.R")
for x,y in zip(foo.x, foo.y):
	plt += expand(expand(plt__sid_sub_, out1=x, plt=y), sid=["231", "232"], sub="epi")

plt__sid_sub__sid_sub = "plts/{out1},{out2},{plt},{sid},{sub}.pdf"
plt__sid_sub__sid_sub_ = "plts/{out1},{out2},{plt},{{sid}},{{sub}}.pdf"
foo = glob_wildcards("code/10-plt__sid_sub__sid_sub-{x},{y},{z}.R")
for x,y,z in zip(foo.x, foo.y, foo.z):
	plt += expand(plt__sid_sub__sid_sub_, out1=x, out2=y, plt=z)

plt__sid__one_sid_all_sub = "plts/{out1},{out2},{plt},{sid}.pdf"
plt__sid__one_sid_all_sub_ = "plts/{out1},{out2},{plt},{{sid}}.pdf"
foo = glob_wildcards("code/10-plt__sid__one_sid_all_sub-{x},{y},{z}.R")
for x,y,z in zip(foo.x, foo.y, foo.z):
	plt += expand(plt__sid__one_sid_all_sub_, out1=x, out2=y, plt=z)

plt__all_sid_all_sub = "plts/{out1},{plt}.pdf"
foo = glob_wildcards("code/10-plt__all_sid_all_sub-{x},{y}.R")
for x,y in zip(foo.x, foo.y):
	plt += expand(plt__all_sid_all_sub, out1=x, plt=y)

plt__all_sid_one_sub = "plts/{out1},{plt},{sub}.pdf"
plt__all_sid_one_sub_ = "plts/{out1},{plt},{{sub}}.pdf"
foo = glob_wildcards("code/10-plt__all_sid_one_sub-{x},{y}.R")
for x,y in zip(foo.x, foo.y):
	plt += expand(plt__all_sid_one_sub_, out1=x, plt=y)

plt__all_sid__all_sid_all_sub = "plts/{out1},{out2},{plt}.pdf"
foo = glob_wildcards("code/10-plt__all_sid__all_sid_all_sub-{x},{y},{z}.R")
for x,y,z in zip(foo.x, foo.y, foo.z):
	plt += expand(plt__all_sid__all_sid_all_sub, out1=x, out2=y, plt=z)

# pat = re.compile(r'^((?!roi).)*$')
# plt = [p for p in plt if pat.match(p)]

pat = re.compile(r'^((?!trj).)*$')
plt = [p for p in plt if pat.match(p)]

# pat = re.compile(r'^.*trj.*$')
# qlt = [p for p in plt if pat.match(p)]
# plt += expand(qlt, sid=SID, sub="epi")
# pat = re.compile(r'^((?!trj.*{).)*$')
# plt = [p for p in plt if pat.match(p)]

# visuals that require so many inputs, 
# they don't fit with the above schema...
qlt = []

# transition crypts
def qlt_tcs_lvx_r(n): return("code/10-qlt_tcs-lv"+str(n)+",{x}.R")
def qlt_tcs_lvx_p(n): return("plts/tcs,lv"+str(n)+",{x}.pdf")
for n in [1, 2]:
    foo = glob_wildcards(qlt_tcs_lvx_r(n))
    for x in foo.x: 
        qlt += expand(qlt_tcs_lvx_p(n), x=x)

rule all:
	input:
		expand([raw, fil, pol, pro, roi, ccc, sig, ist, lv1], sid=SID),
		expand([sub, jst, lv2], sid=sid, sub=SUB),
		#expand([rep, trj], sid=SID, sub=["epi"]),
		expand([pbs], sid=SID, sub=SUB),
		expand(plt, sid=SID, sub=SUB),
		expand([ctx], sid=SID), qlt

# analysis =========================================

# coercion
for did in DID:
	rule:
		priority: 99
		input:	"code/01-raw.R"
		params:	expand("data/raw/{did}", did=did)[0], md, ex
		output:	directory(expand(raw, sid=MD["sid"][MD["did"] == did]))
		log:    expand("logs/raw-{did}.Rout", did=did)
		shell: '''R CMD BATCH\\
		--no-restore --no-save "--args wcs={wildcards}\
		{params} {output}" {input[0]} {log}'''

# filtering
def foo(wildcards): return(expand(
	"data/raw/{did}/positions.csv.gz", 
	did=MD["did"][MD["sid"] == wildcards.sid]))
rule fil:
	priority: 98
	input:	"code/02-fil.R", raw, foo
	output:	fil
	log:    "logs/fil-{sid}.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {input[2]} {output}" {input[0]} {log}'''

# polygons
def foo(wildcards): return(expand(
	"data/raw/{did}/polygons.csv.gz", 
	did=MD["did"][MD["sid"] == wildcards.sid]))
rule pol:
	input: 	"code/03-pol.R", fil, foo
	output:	pol
	log: 	"logs/pol-{sid}.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {input[2]} {output}" {input[0]} {log}'''

# processing
rule pro:
	threads: 20
	priority: 96
	input:	"code/03-pro.R", fil, pbs_lv1
	output:	pro
	log:    "logs/pro-{sid}.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {input[2]} {output} ths={threads}" {input[0]} {log}'''

# regions
def foo(wildcards): return(expand(
	"imgs/{did}/shapes/{sid}/{roi}",
	did=MD["did"][MD["sid"] == wildcards.sid],
	sid=wildcards.sid, roi=ROI[wildcards.sid]))
rule roi:
	input:	"code/03-roi.R", fil, x = foo
	params: lambda wc, input: ";".join(input.x)
	output:	roi
	log:	"logs/roi-{sid}.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {params} {output[0]}" {input[0]} {log}'''

# communication
rule ccc:
	threads: 40
	input:	"code/04-ccc.R", fil
	output:	ccc
	log:    "logs/ccc-{sid}.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {output[0]} ths={threads}" {input[0]} {log}'''

# signatures
rule sig:
	priority: 97
	threads: 10
	input:	"code/03-sig.R", fil
	output:	sig
	log:    "logs/sig-{sid}.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {output[0]} ths={threads}" {input[0]} {log}'''

# clustering
rule ist:
	priority: 97
	input:	"code/03-ist.R", fil, pbs_lv1
	output:	ist
	log:    "logs/ist-{sid}.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {input[2]} {output}" {input[0]} {log}'''	

# labelling
rule lv1:
	priority: 98
	input:	"code/00-lab.R", ist,
			"meta/lab/lv1.json"
	output:	lv1
	log:    "logs/lv1-{sid}.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {input[2]} {output}" {input[0]} {log}'''	

# subsetting
rule sub:
	priority: 97
	input:	"code/04-sub.R", fil, lv1
	output:	expand("outs/sub-{{sid}},{sub}.rds", sub=SUB)
	log:    "logs/sub-{sid}.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {input[2]} {output}" {input[0]} {log}'''	

# subclustering
rule jst:
	priority: 96
	input:	"code/05-jst.R", sub, pbs_lv2
	output:	jst
	log:    "logs/jst-{sid},{sub}.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {input[2]} {output}" {input[0]} {log}'''	

# labelling
rule lv2:
	priority: 95
	input:	"code/00-lab.R", jst,
			"meta/lab/lv2,{sub}.json"
	output:	lv2
	log:    "logs/lv2-{sid},{sub}.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {input[2]} {output}" {input[0]} {log}'''

# contexts
rule ctx:
	priority: 95
	input:	"code/06-ctx.R", 
			x = expand(fil, sid=SID),
			y = expand(jst, sid=SID, sub=SUB)
	params:	lambda wc, input: ";".join(input.x),
			lambda wc, input: ";".join(input.y)
	output:	expand(ctx, sid=SID)
	log:    "logs/ctx.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{params[0]} {params[1]} {output}" {input[0]} {log}'''	

# reprocessing
rule rep:
	wildcard_constraints: sub = "epi"
	threads: 20
	priority: 94
	input:	"code/05-rep.R", sub, pbs_lv2
	output:	rep
	log:    "logs/rep-{sid},{sub}.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {input[2]} {output} ths={threads}" {input[0]} {log}'''

# trajectory
rule trj:
	priority: 92
	wildcard_constraints: sub = "epi"
	input:	"code/06-trj.R", rep, jst
	output:	trj
	log:    "logs/trj-{sid},{sub}.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {input[2]} {output}" {input[0]} {log}'''	

# # intergation
# rule add:
# 	wildcard_constraints: sub = "epi"
# 	threads: 80
# 	priority: 93
# 	input:	"code/06-red.R",
# 			x = expand("outs/rep-{sid},{{sub}}.rds", sid=SID)
# 	params: lambda wc, input: ";".join(input.x)
# 	output:	directory(add)
# 	log:    "logs/add-{sub}.Rout"
# 	shell: '''R CMD BATCH\\
# 	--no-restore --no-save "--args wcs={wildcards}\
# 	{params} {output} ths={threads}" {input[0]} {log}'''	

# # trajectory
# rule trj:
# 	priority: 92
# 	wildcard_constraints: sub = "epi"
# 	input:	"code/07-trj.R", add,
# 			x = expand("outs/jst-{sid},{{sub}}.rds", sid=SID)
# 	params: lambda wc, input: ";".join(input.x)
# 	output:	trj
# 	log:    "logs/trj-{sub}.Rout"
# 	shell: '''R CMD BATCH\\
# 	--no-restore --no-save "--args wcs={wildcards}\
# 	{input[1]} {params} {output}" {input[0]} {log}'''	

# profiles	
rule pbs:
	priority: 94
	threads: 20
	input:	"code/00-pbs.R", sub, lv2
	output:	pbs
	log:    "logs/pbs-{sid},{sub}.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {input[2]} {output} ths={threads}" {input[0]} {log}'''	

# # pooling
# rule add:
# 	priority: 96
# 	input:	"code/05-add.R",
# 			x = expand("outs/sub-{sid},{{sub}}.rds", sid=SID)
# 	params: lambda wc, input: ";".join(input.x)
# 	output:	directory(add)
# 	log:    "logs/add-{sub}.Rout"
# 	shell: '''R CMD BATCH\\
# 	--no-restore --no-save "--args wcs={wildcards}\
# 	{params} {output}" {input[0]} {log}'''	

# # reprocessing
# rule rep:
# 	threads: 30
# 	priority: 95
# 	input:	"code/06-rep.R", add
# 	output:	rep
# 	log:    "logs/rep-{sub}.Rout"
# 	shell: '''R CMD BATCH\\
# 	--no-restore --no-save "--args wcs={wildcards}\
# 	{input[1]} {output} ths={threads}" {input[0]} {log}'''

# # integration
# rule red:
# 	threads: 30
# 	priority: 94
# 	input:	"code/07-red.R", rep
# 	output:	red
# 	log:    "logs/red-{sub}.Rout"
# 	shell: '''R CMD BATCH\\
# 	--no-restore --no-save "--args wcs={wildcards}\
# 	{input[1]} {output} ths={threads}" {input[0]} {log}'''

# # subclustering
# rule clu:
# 	threads: 30
# 	priority: 93
# 	input:	"code/08-clu.R", red
# 	params: lambda wc: {"imm":0.6, "epi":1.2, "str":0.8}[wc.sub]
# 	output:	clu
# 	log:    "logs/clu-{sub}.Rout"
# 	shell: '''R CMD BATCH\\
# 	--no-restore --no-save "--args wcs={wildcards}\
# 	{input[1]} {output} ths={threads} res={params}" {input[0]} {log}'''	

# pooling
rule add_by_sid:
	priority: 90
	input:	"code/05-add.R", fil, ist, lv1, 
			x = expand("outs/jst-{{sid}},{sub}.rds", sub=SUB)
	params:	lambda wc, input: ";".join(input.x)
	output:	add_by_sid
	log:    "logs/add_by_sid-{sid}.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {input[2]} {input[3]} {params}\
	{output}" {input[0]} {log}'''	

rule add_by_sub:
	priority: 90
	input:	"code/05-add.R", fil, ist, lv1, 
			x = expand("outs/jst-{sid},{{sub}}.rds", sid=SID)
	params:	lambda wc, input: ";".join(input.x)
	output:	add_by_sub
	log:    "logs/add_by_sub-{sub}.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {input[2]} {input[3]} {params}\
	{output}" {input[0]} {log}'''	

# plotting =========================================

rule plt_raw:
	priority: 9
	input:	"code/10-plt_raw,{out2}-{plt}.R",
			"outs/raw-{sid}/se.rds", 
			"outs/{out2}-{sid}.rds"
	output:	plt_raw
	log:	"logs/plt-raw,{out2}-{plt},{sid}.Rout"
	shell: '''R CMD BATCH\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {input[2]} {output}" {input[0]} {log}'''

a_sid = "outs/{out1}-{sid}.rds"
b_sid = "outs/{out2}-{sid}.rds"
a_sid_sub = "outs/{out1}-{sid},{sub}.rds"
b_sid_sub = "outs/{out2}-{sid},{sub}.rds"
a_all_sid = expand("outs/{{out1}}-{sid}.rds", sid=SID)
b_all_sid = expand("outs/{{out2}}-{sid}.rds", sid=SID)
a_all_sid_one_sub = expand("outs/{{out1}}-{sid},{{sub}}.rds", sid=SID)
b_all_sid_one_sub = expand("outs/{{out2}}-{sid},{{sub}}.rds", sid=SID)
a_one_sid_all_sub = expand("outs/{{out1}}-{{sid}},{sub}.rds", sub=SUB)
b_one_sid_all_sub = expand("outs/{{out2}}-{{sid}},{sub}.rds", sub=SUB)
a_all_sid_all_sub = expand("outs/{{out1}}-{sid},{sub}.rds", sid=SID, sub=SUB)
b_all_sid_all_sub = expand("outs/{{out2}}-{sid},{sub}.rds", sid=SID, sub=SUB)

rule plt__all_sid:
	priority: 8
	input:	"code/10-plt__all_sid-{out1},{plt}.R", x = a_all_sid
	log:	"logs/plt__all_sid-{out1},{plt}.Rout"
	params:	lambda wc, input: ";".join(input.x)
	output:	plt__all_sid
	shell: '''R CMD BATCH\
	--no-restore --no-save "--args wcs={wildcards}\
	{params} {output}" {input[0]} {log}'''

rule plt__sid:
	priority: 8
	input:	"code/10-plt__sid-{out1},{plt}.R", a_sid
	log:	"logs/plt__sid-{out1},{plt},{sid}.Rout"
	output:	plt__sid
	shell: '''R CMD BATCH\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {output}" {input[0]} {log}'''

rule plt__sid__sid:
	priority: 8
	input:	"code/10-plt__sid__sid-{out1},{out2},{plt}.R", a_sid, b_sid
	log:	"logs/plt__sid__sid-{out1},{out2},{plt},{sid}.Rout"
	output:	plt__sid__sid
	shell: '''R CMD BATCH\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {input[2]} {output}" {input[0]} {log}'''

rule plt__all_sid__all_sid:
	priority: 8
	input:	"code/10-plt__all_sid__all_sid-{out1},{out2},{plt}.R", x = a_all_sid, y = b_all_sid
	log:	"logs/plt__all_sid__all_sid-{out1},{out2},{plt}.Rout"
	params:	lambda wc, input: ";".join(input.x),
			lambda wc, input: ";".join(input.y)
	output:	plt__all_sid__all_sid
	shell: '''R CMD BATCH\
	--no-restore --no-save "--args wcs={wildcards}\
	{params} {output}" {input[0]} {log}'''

rule plt__all_sid_all_sub:
	priority: 8
	input:	"code/10-plt__all_sid_all_sub-{out1},{plt}.R", x = a_all_sid_all_sub
	log:	"logs/plt__all_sid_all_sub-{out1},{plt}.Rout"
	params:	lambda wc, input: ";".join(input.x)
	output:	plt__all_sid_all_sub
	shell: '''R CMD BATCH\
	--no-restore --no-save "--args wcs={wildcards}\
	{params} {output}" {input[0]} {log}'''

rule plt__all_sid_one_sub:
	priority: 8
	input:	"code/10-plt__all_sid_one_sub-{out1},{plt}.R", x = a_all_sid_one_sub
	log:	"logs/plt__all_sid_one_sub-{out1},{plt},{sub}.Rout"
	params:	lambda wc, input: ";".join(input.x)
	output:	plt__all_sid_one_sub
	shell: '''R CMD BATCH\
	--no-restore --no-save "--args wcs={wildcards}\
	{params} {output}" {input[0]} {log}'''

rule plt__sid_sub:
	priority: 8
	input:	"code/10-plt__sid_sub-{out1},{plt}.R", a_sid_sub
	log:	"logs/plt__sid_sub-{out1},{plt},{sid},{sub}.Rout"
	output:	plt__sid_sub
	shell: '''R CMD BATCH\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {output}" {input[0]} {log}'''

rule plt__sid_sub__sid_sub:
	priority: 8
	input:	"code/10-plt__sid_sub__sid_sub-{out1},{out2},{plt}.R", a_sid_sub, b_sid_sub
	log:	"logs/plt__sid_sub__sid_sub-{out1},{out2},{plt},{sid},{sub}.Rout"
	output:	plt__sid_sub__sid_sub
	shell: '''R CMD BATCH\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {input[2]} {output}" {input[0]} {log}'''

rule plt__sid__one_sid_all_sub:
	priority: 8
	input:	"code/10-plt__sid__one_sid_all_sub-{out1},{out2},{plt}.R", a_sid, x = b_one_sid_all_sub
	log:	"logs/plt_all_sub_one_sid-{out1},{out2},{plt},{sid}.Rout"
	params:	lambda wc, input: ";".join(input.x)
	output:	plt__sid__one_sid_all_sub
	shell: '''R CMD BATCH\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {params} {output}" {input[0]} {log}'''

rule plt__all_sid__all_sid_all_sub:
	priority: 8
	input:	"code/10-plt__all_sid__all_sid_all_sub-{out1},{out2},{plt}.R", x = a_all_sid, y = b_all_sid_all_sub
	log:	"logs/plt__all_sid__all_sid_all_sub-{out1},{out2},{plt}.Rout"
	params:	lambda wc, input: ";".join(input.x),
			lambda wc, input: ";".join(input.y)
	output:	plt__all_sid__all_sid_all_sub
	shell: '''R CMD BATCH\
	--no-restore --no-save "--args wcs={wildcards}\
	{params} {output}" {input[0]} {log}'''

# visuals that require so many inputs, 
# they don't fit with the above schema...

def out_sid(out, typ="rds"): return(expand("outs/{out}-{sid}.{typ}", 
	out=out, sid=SID, typ="parquet" if out == "pol" else "rds"))
def out_sid_sub(out): return(expand("outs/{out}-{sid},{sub}.rds", out=out, sid=SID, sub=SUB))
def out_ist(n): return(out_sid("lv1") if n == 1 else out_sid_sub("lv2"))

# transition crypts
for n in [1, 2]:
    rule:
        input:	qlt_tcs_lvx_r(n),
                a = out_sid("fil"),
                b = out_sid("roi"),
                c = out_sid("pol", "parquet"),
                d = out_ist(n)
        log:	"logs/qlt_tcs-lv"+str(n)+",{x}.Rout"
        params:	lambda wc, input: ";".join(input.a),
                lambda wc, input: ";".join(input.b),
                lambda wc, input: ";".join(input.c),
                lambda wc, input: ";".join(input.d)
        output:	qlt_tcs_lvx_p(n)
        shell: '''R CMD BATCH\
        --no-restore --no-save "--args\
        {params} {output}" {input[0]} {log}'''
