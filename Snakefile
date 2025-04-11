# WTx-CosMx SMI on TVA/CRC
import pandas as pd
import re

ex = "meta/exclude.rds" # barcodes overlapping with CARD8 barcode in 1st run
md = "data/raw/metadata.txt" # slide metadata (including FOV-to-section mapping)
MD = pd.read_csv(md, dtype="string", delimiter="\t")

#onstart: shell("Rscript meta/sig/sig.R")

DID = []
for i in range(len(MD)):
	did = MD["did"][i]
	if did not in DID:
		DID += [did]
SID = sorted(set(MD["sid"]))
SID = [x for x in SID if x != "241"]
SUB = ["epi", "imm", "str"]

# regions of interest
ROI = {}
rmv = ["110_ROI4_RC", "120_ROI7_RC", "120_ROI8_RC", 
    "120_ROI9_RC", "221_ROI1_RT2", "231_ROI11_TC"]
# these exhibit transitions in H&E,
# but not the adjacent section's IF
for sid in SID+["241"]:
	did = MD["did"][MD["sid"] == sid]
	foo = "imgs/{did}/shapes/{sid}/{{roi}}.pickle"
	foo = expand(foo, did=did, sid=sid)[0]
	foo = glob_wildcards(foo).roi
	foo = [x for x in foo if x not in rmv]
	if sid not in ROI.keys():
		ROI[sid] = {}
	ROI[sid] = foo

# snPATHO-seq reference profiles
pbs_lv1 = "data/ref/outs/pbs-lv1.rds"
pbs_lv2 = "data/ref/outs/pbs-lv2,{sub}.rds"

# targets ==========================================

# processing
raw = "outs/raw-{sid}"
fil = "outs/fil-{sid}.rds"
roi = "outs/roi-{sid}.rds"
pol = "outs/pol-{sid}.parquet"
# downstream
pro = "outs/pro-{sid}.rds"
clu = "outs/clu-{sid}.rds"
sig = "outs/sig-{sid}.rds"
ccc = "outs/ccc-{sid}.rds"
# clustering
ist = "outs/ist-{sid}.rds"
lv1 = "outs/lv1-{sid}.rds"
pbs = "outs/pbs.rds"
# subclustering
sub = "outs/sub-{sid},{sub}.rds"
jst = "outs/jst-{sid},{sub}.rds"
lv2 = "outs/lv2-{sid},{sub}.rds"
# epithelia
rep = "outs/rep-{sid}.rds"
trj = "outs/trj-{sid}.rds"
kst = "outs/kst-{sid}.rds"
ctx = "outs/ctx-{sid}.rds"
qbs = "outs/qbs-{sub}.rds"
pbs_kst = "outs/pbs_kst.rds"
mgs = "outs/mgs-{sub}.rds"

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

# one by nan
pdf = "plts/{out},{plt},{ids}.pdf"
foo = "code/10-plt__{x}-{{a}},{{p}}.R"
for x in ["sid", "sub"]:
    i = {"sid": SID, "sub": SUB}[x]
    bar = glob_wildcards(expand(foo, x=x)[0])
    for a,p in zip(bar.a, bar.p):
        plt += expand(pdf, out=a, plt=p, ids=i)

# all by nan
pdf = "plts/{out},{plt}.pdf"
foo = "code/10-plt__{x}-{{a}},{{p}}.R"
for x in ["all_sid", "all_sub", "all_sid_all_sub"]:
    bar = glob_wildcards(expand(foo, x=x)[0])
    for a,p in zip(bar.a, bar.p):
        plt += expand(pdf, out=a, plt=p)

# two all by nan
pdf = "plts/{out1},{out2},{plt}.pdf"
foo = "code/10-plt__{x}__{y}-{{a}},{{b}},{{p}}.R"
for xy in [
	["all_sid", "all_sid"],
	["all_sid", "all_sid_all_sub"]]:
    bar = glob_wildcards(expand(foo, x=xy[0], y=xy[1])[0])
    for a,b,p in zip(bar.a, bar.b, bar.p):
        plt += expand(pdf, out1=a, out2=b, plt=p)

# two by sid
pdf = "plts/{out1},{out2},{plt},{{sid}}.pdf"
foo = "code/10-plt__{x}__{y}-{{a}},{{b}},{{p}}.R"
for xy in [
    ["sid", "sid"],
    ["sid", "one_sid_all_sub"],
    ["one_sid_all_sub", "one_sid_all_sub"]]:
    bar = glob_wildcards(expand(foo, x=xy[0], y=xy[1])[0])
    for a,b,p in zip(bar.a, bar.b, bar.p):
        plt += expand(pdf, out1=a, out2=b, plt=p)

# one by sub
pdf = "plts/{out1},{plt},{{sub}}.pdf"
foo = glob_wildcards("code/10-plt__all_sid_one_sub-{x},{y}.R")
for x,y in zip(foo.x, foo.y):
	plt += expand(expand(pdf, out1=x, plt=y), sub=SUB)

# one by sid-sub
pdf = "plts/{out1},{plt},{{sid}},{{sub}}.pdf"
foo = glob_wildcards("code/10-plt__sid_sub-{x},{y}.R")
for x,y in zip(foo.x, foo.y):
	plt += expand(expand(pdf, out1=x, plt=y), sid=SID, sub="epi")

# two by sid-sub
pdf = "plts/{out1},{out2},{plt},{{sid}},{{sub}}.pdf"
foo = glob_wildcards("code/10-plt__sid_sub__sid_sub-{x},{y},{z}.R")
for x,y,z in zip(foo.x, foo.y, foo.z):
	plt += expand(pdf, out1=x, out2=y, plt=z)

#pat = re.compile(r'^((?!trj).)*$')
#plt = [p for p in plt if pat.match(p)]
#qlt = [q for q in qlt if pat.match(q)]

pat = re.compile(r'^.*kst.*$')
qlt = [p for p in plt if pat.match(p)]
plt += expand(qlt, sid=SID, sub="epi")

# visuals that require so many inputs, 
# they don't fit with the above schema...
qlt = []

# epithelial trajectors
qlt_trj_r = "code/10-qlt_trj-{x},{p}.R"
qlt_trj_p = "plts/trj,{x},{p},{{sid}}.pdf"
foo = glob_wildcards(qlt_trj_r)
for x,p in zip(foo.x, foo.p):
    qlt += expand(qlt_trj_p, x=x, p=p)

# transition crypts
qlt_tcs_r = "code/10-qlt_tcs-{x},{p}.R"
qlt_tcs_p = "plts/tcs,{x},{p}.pdf"
foo = glob_wildcards(qlt_tcs_r)
for x,p in zip(foo.x, foo.p):
    qlt += expand(qlt_tcs_p, x=x, p=p)

rule all:
    input:
        expand([raw, fil, pol, roi], sid=SID+["241"]),
        expand([clu, pro, ccc, sig], sid=SID),
        expand([ist, lv1, pbs, ctx], sid=SID),
        expand([sub, jst, lv2, qbs, mgs], sid=sid, sub=SUB),
        expand([kst, rep, trj], sid=SID, sub="epi"),
        expand(plt, sid=SID, sub=SUB),
        expand(qlt, sid=SID)#, pbs_kst

# analysis =========================================

# reference
rule ref:
    priority: 99
    input:  "code/00-ref.R", "data/gca.rds"
    output: "outs/gca.rds"
    log:    "logs/ref.Rout"
    shell: '''R CMD BATCH\\
    --no-restore --no-save "--args wcs={wildcards}\
    {input[1]} {output}" {input[0]} {log}'''

# CNV clustering
rule clu:
	priority: 99
	threads: 10
	input:	"code/01-clu.R", "outs/cnv-{sid}.rds"
	output:	clu
	log:    "logs/clu-{sid}.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards} 
	{input[1]} {output} ths={threads}" {input[0]} {log}'''
	
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
	
# signatures
rule sig:
    priority: 97
    threads: 10
    input:  "code/03-sig.R", fil, 
            x = expand("meta/sig/sig-{sub}.json", sub=SUB)
    params:	lambda wc, input: ";".join(input.x)
    output:	sig
    log:    "logs/sig-{sid}.Rout"
    shell: '''R CMD BATCH\\
    --no-restore --no-save "--args wcs={wildcards}\
    {input[1]} {params} {output} ths={threads}" {input[0]} {log}'''

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

# reclustering
rule kst:
    priority: 95
    input:	"code/06-kst.R", rules.ref.output, 
            "outs/sub-{sid},epi.rds", 
            "outs/jst-{sid},epi.rds"
    output:	kst
    log:    "logs/kst-{sid}.Rout"
    shell: '''R CMD BATCH\\
    --no-restore --no-save "--args wcs={wildcards}\
    {input[1]} {input[2]} {input[3]} {output}" {input[0]} {log}'''	

# profiles	
one = expand("outs/fil-{sid}.rds", sid=SID)
two = expand("outs/ist-{sid}.rds", sid=SID)
rule pbs:
	priority: 94
	threads: 20
	input:	"code/00-pbs.R", x=one, y=two
	params:	lambda wc, input: ";".join(input.x),
			lambda wc, input: ";".join(input.y)
	output:	pbs
	log:    "logs/pbs.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{params} {output} ths={threads}" {input[0]} {log}'''	

one = expand("outs/sub-{sid},{{sub}}.rds", sid=SID)
two = expand("outs/lv2-{sid},{{sub}}.rds", sid=SID)
rule qbs:
	threads: 30
	priority: 95
	input:	"code/00-pbs.R", x=one, y=two
	params:	lambda wc, input: ";".join(input.x),
			lambda wc, input: ";".join(input.y)
	output:	qbs
	log:    "logs/qbs-{sub}.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{params} {output} ths={threads}" {input[0]} {log}'''

rule pbs_kst:
	threads: 30
	priority: 95
	input:	"code/00-pbs.R", 
            x=expand("outs/sub-{sid},epi.rds", sid=SID), 
            y=expand("outs/kst-{sid}.rds", sid=SID)
	params:	lambda wc, input: ";".join(input.x),
			lambda wc, input: ";".join(input.y)
	output:	pbs_kst
	log:    "logs/mgs_kst.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{params} {output} ths={threads}" {input[0]} {log}'''

one = expand("outs/sub-{sid},{{sub}}.rds", sid=SID)
two = expand("outs/lv2-{sid},{{sub}}.rds", sid=SID)
rule mgs:
	threads: 30
	priority: 95
	input:	"code/00-mgs.R", x=one, y=two
	params:	lambda wc, input: ";".join(input.x),
			lambda wc, input: ";".join(input.y)
	output:	mgs
	log:    "logs/mgs-{sub}.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{params} {output} ths={threads}" {input[0]} {log}'''

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
	threads: 20
	priority: 94
	input:	"code/05-rep.R", 
            "outs/sub-{sid},epi.rds",
            expand(pbs_lv2, sub="epi")
	output:	rep
	log:    "logs/rep-{sid}.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {input[2]} {output} ths={threads}" {input[0]} {log}'''

# trajectory
rule trj:
	priority: 92
	input:	"code/06-trj.R", rep, roi, "outs/jst-{sid},epi.rds"
	output:	trj
	log:    "logs/trj-{sid}.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {input[2]} {input[3]} {output}" {input[0]} {log}'''	

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
a_sub = "outs/{out1}-{sub}.rds"
b_sub = "outs/{out2}-{sub}.rds"
a_sid_sub = "outs/{out1}-{sid},{sub}.rds"
b_sid_sub = "outs/{out2}-{sid},{sub}.rds"
a_all_sub = expand("outs/{{out1}}-{sub}.rds", sub=SUB)
b_all_sub = expand("outs/{{out2}}-{sub}.rds", sub=SUB)
a_all_sid = expand("outs/{{out1}}-{sid}.rds", sid=SID)
b_all_sid = expand("outs/{{out2}}-{sid}.rds", sid=SID)
a_all_sid_one_sub = expand("outs/{{out1}}-{sid},{{sub}}.rds", sid=SID)
b_all_sid_one_sub = expand("outs/{{out2}}-{sid},{{sub}}.rds", sid=SID)
a_one_sid_all_sub = expand("outs/{{out1}}-{{sid}},{sub}.rds", sub=SUB)
b_one_sid_all_sub = expand("outs/{{out2}}-{{sid}},{sub}.rds", sub=SUB)
a_all_sid_all_sub = expand("outs/{{out1}}-{sid},{sub}.rds", sid=SID, sub=SUB)
b_all_sid_all_sub = expand("outs/{{out2}}-{sid},{sub}.rds", sid=SID, sub=SUB)	

def out(by, n=None):
    o = "out" if n is None else "out"+str(n)
    d = {
        "sid": "outs/{"+o+"}-{x}.rds",
        "sub": "outs/{"+o+"}-{x}.rds",
        "sid_sub": "outs/{"+o+"}-{x}.rds",
        "all_sid": expand("outs/{{"+o+"}}-{x}.rds", x=SID),
        "all_sub": expand("outs/{{"+o+"}}-{x}.rds", x=SUB),
        "all_sid_one_sub": expand("outs/{{"+o+"}}-{y},{{x}}.rds", y=SID),
		"one_sid_all_sub": expand("outs/{{"+o+"}}-{{x}},{y}.rds", y=SUB),
        "all_sid_all_sub": expand("outs/{{"+o+"}}-{x},{y}.rds", x=SID, y=SUB)
    }
    return(d[by])

def collect(x): return(x if type(x) == str else ";".join(x))

pat = "plt__{by}-{{out}},{{plt}}"
for by in ["sid", "sub", "all_sid_one_sub"]:
    rule:
        name:   "plt__%s" % by
        input:  expand("code/10-"+pat+".R", by=by), a = out(by)
        params: lambda wc, input: collect(input.a)
        log:    expand("logs/"+pat+",{{x}}.Rout", by=by) 
        output: "plts/{out},{plt},{x}.pdf"
        shell: '''R CMD BATCH\
        --no-restore --no-save "--args wcs={wildcards}\
        {params} {output}" {input[0]} {log}'''

pat = "plt__{by}-{{out}},{{plt}}"
for by in ["all_sid", "all_sub", "all_sid_all_sub"]:
    rule:
        name:   "plt__%s" % by
        input:	expand("code/10-"+pat+".R", by=by), a = out(by)
        params: lambda wc, input: collect(input.a)
        log:	expand("logs/"+pat+".Rout", by=by)
        output:	"plts/{out},{plt}.pdf"
        shell: '''R CMD BATCH\
        --no-restore --no-save "--args wcs={wildcards}\
        {params} {output}" {input[0]} {log}'''

pat = "plt__{by1}__{by2}-{{out1}},{{out2}},{{plt}}"
for by in [
	["all_sid","all_sid"],
	["all_sid","all_sid_all_sub"]]:
    rule:
        name:   "plt__%s__%s" % (by[0], by[1])
        input:  expand("code/10-"+pat+".R", by1=by[0], by2=by[1]),
                a = out(by[0], 1), b = out(by[1], 2)
        params: lambda wc, input: ";".join(input.a),
                lambda wc, input: ";".join(input.b)
        log:    expand("logs/"+pat+".Rout", by1=by[0], by2=by[1])
        output:	"plts/{out1},{out2},{plt}.pdf"
        shell: '''R CMD BATCH\
        --no-restore --no-save "--args wcs={wildcards}\
        {params} {output}" {input[0]} {log}'''

# plt by sid
pat = "plt__{by1}__{by2}-{{out1}},{{out2}},{{plt}}"
for by in [
    ["sid","sid"],
    ["sid","one_sid_all_sub"],
    ["one_sid_all_sub","one_sid_all_sub"]]:
    rule:
        name:   "plt__%s__%s" % (by[0], by[1])
        input:  expand("code/10-"+pat+".R", by1=by[0], by2=by[1]),
                a = out(by[0], 1), b = out(by[1], 2)
        params: lambda wc, input: collect(input.a),
                lambda wc, input: collect(input.b)
        log:    expand("logs/"+pat+",{{x}}.Rout", by1=by[0], by2=by[1])
        output: "plts/{out1},{out2},{plt},{x}.pdf"
        shell: '''R CMD BATCH\
        --no-restore --no-save "--args wcs={wildcards}\
        {params} {output}" {input[0]} {log}'''

pat = "plt__{by1}__{by2}-{{out1}},{{out2}},{{plt}}"
for by in [
	["sid_sub","sid_sub"]]:
    rule:
        name:   "plt__%s__%s" % (by[0], by[1])
        input:  expand("code/10-"+pat+".R", by1=by[0], by2=by[1]),
                a = out(by[0], 1), b = out(by[1], 1)
        params: lambda wc, input: collect(input.a),
                lambda wc, input: collect(input.b)
        log:    expand("logs/"+pat+",{{x}},{{y}}.Rout", by1=by[0], by2=by[1])
        output: "plts/{out1},{out2},{plt},{x},{y}.pdf"
        shell: '''R CMD BATCH\
        --no-restore --no-save "--args wcs={wildcards}\
        {params} {output}" {input[0]} {log}'''

# visuals that require so many inputs, 
# they don't fit with the above schema

def out_sid(out, typ="rds"): return(expand("outs/{out}-{sid}.{typ}", 
	out=out, sid=SID, typ="parquet" if out == "pol" else "rds"))
def out_sid_sub(out): return(expand("outs/{out}-{sid},{sub}.rds", out=out, sid=SID, sub=SUB))
def out_foo(x):
    fun = out_sid
    if x in ["jst", "lv2"]:
        fun = out_sid_sub
    return(fun(x))

# epithelial trajectirs
def xxx(x):
	if x in ["sub", "jst", "lv2"]:
		rds = "outs/{x}-{{sid}},{sub}.rds"
	else:
		rds = "outs/{x}-{{sid}}.rds"
	return(expand(rds, x=x, sub=SUB))

qlt_trj_r = "code/10-qlt_trj-{x},{p}.R"
qlt_trj_p = "plts/trj,{x},{p},{{sid}}.pdf"

foo = glob_wildcards(qlt_trj_r)
for x,p in zip(foo.x, foo.p):
    rule:
        priority: 7
        input:  expand(qlt_trj_r, x=x, p=p), trj, fil, z = xxx(x)
        log:    "logs/qlt_trj-"+x+","+p+",{sid}.Rout"
        params:	lambda wc, input: ";".join(input.z)
        output:	expand(qlt_trj_p, x=x, p=p)
        shell: '''R CMD BATCH\
        --no-restore --no-save "--args wcs={wildcards}\
        {input[1]} {input[2]} {params} {output}" {input[0]} {log}'''

# transition crypts
foo = glob_wildcards(qlt_tcs_r)
for x,p in zip(foo.x, foo.p):
    rule:
        priority: 7
        input:	expand(qlt_tcs_r, x=x, p=p),
                a = out_sid("fil"),
                b = out_sid("roi"),
                c = out_sid("pol", "parquet"),
                d = out_foo(x)
        log:	"logs/qlt_tcs-"+x+","+p+".Rout"
        params:	lambda wc, input: ";".join(input.a),
                lambda wc, input: ";".join(input.b),
                lambda wc, input: ";".join(input.c),
                lambda wc, input: ";".join(input.d)
        output:	expand(qlt_tcs_p, x=x, p=p)
        shell: '''R CMD BATCH\
        --no-restore --no-save "--args\
        {params} {output}" {input[0]} {log}'''
