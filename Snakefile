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
SIP = sorted(set(MD["sid"]))
SID = [x for x in SIP if x != "241"]
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
ctx = "outs/ctx-{sid}.rds"
cty = "outs/cty-{sid}.rds"
ccc = "outs/ccc-{sid}.rds"
pro = "outs/pro-{sid}.rds"
sig = "outs/sig-{sid}.rds"
clu = "outs/clu-{sid}.rds"
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
req = "outs/req-{sid}.rds"
mgs = "outs/mgs-{sub}.rds"
qbs = "outs/qbs-{sub}.rds"
qbt = "outs/qbt.rds"

# plotting
plt = []

# one
pdf = "plts/{out},{plt}.pdf"
foo = "code/10-plt__one-{a},{p}.R"
bar = glob_wildcards(foo)
for a,p in zip(bar.a, bar.p):
    plt += expand(pdf, out=a, plt=p)

# raw
plt_raw = "plts/raw/raw,{out2},{plt},{sid}.pdf"
plt_raw_ = "plts/raw/raw,{out2},{plt},{{sid}}.pdf"
bar = glob_wildcards("code/10-plt_raw,{y}-{z}.R")
for y,z in zip(bar.y, bar.z):
    plt += expand(plt_raw_, out2=y, plt=z)

# one by nan
pdf = "plts/{out},{plt},{ids}.pdf"
foo = "code/10-plt__{x}-{{a}},{{p}}.R"
for x in ["sip", "sid", "sub"]:
    i = {"sip": SIP, "sid": SID, "sub": SUB}[x]
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
    ["all_raw", "all_sid"],
	["all_sid", "all_sid"],
    ["all_sip", "all_sip"],
	["all_sid", "all_sid_all_sub"]]:
    bar = glob_wildcards(expand(foo, x=xy[0], y=xy[1])[0])
    for a,b,p in zip(bar.a, bar.b, bar.p):
        plt += expand(pdf, out1=a, out2=b, plt=p)

# two by sid
pdf = "plts/{out1},{out2},{plt},{sid}.pdf"
foo = "code/10-plt__{x}__{y}-{{a}},{{b}},{{p}}.R"
for xy in [
    ["sid", "sid"],
    ["sip", "sip"],
    ["sid", "one_sid_all_sub"],
    ["one_sid_all_sub", "one_sid_all_sub"]]:
    bar = glob_wildcards(expand(foo, x=xy[0], y=xy[1])[0])
    i = [SIP if "sip" in xy else SID][0]
    for a,b,p in zip(bar.a, bar.b, bar.p):
        plt += expand(pdf, out1=a, out2=b, plt=p, sid=i)

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

# visuals that require so many inputs, 
# they don't fit with the above schema...
qlt = []

# epithelial trajectors
qlt_trj_r = "code/10-qlt_trj-{x},{p}.R"
qlt_trj_p = "plts/trj,{x},{p}.pdf"
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
        # setup
        expand([raw, fil, pol], sid=SIP),
        expand([roi, sig, lv1], sid=SIP),
        expand([lv2], sid="241", sub="epi"),
        # clustering
        expand([pro, ist, pbs], sid=SID),
        # subclustering
        expand([sub, jst, lv2, qbs, mgs], sid=SID, sub=SUB),
        expand([kst, qbt, req], sid=SID),
        # downstream
        expand([ctx, cty, ccc, rep, trj], sid=SID),
        # visualization
        expand(plt + qlt, sid=SID),
        # collection
        expand("outs/res-{sid}", sid=SID)

# write session info to .txt file
rule inf:
    input:	"code/09-inf.R"
    output:	"inf.txt"
    log:	"logs/inf.Rout" 
    shell:	'''{R} CMD BATCH\\
    --no-restore --no-saves "--args\\
    {output}" {input} {log}'''

# analysis =========================================

def pool(x): return(x if type(x) == str else ";".join(x))

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

# polygons
def foo(wildcards): return(expand(
	"data/raw/{did}/polygons.csv.gz", 
	did=MD["did"][MD["sid"] == wildcards.sid]))
rule pol:
    priority: 97
    input: 	"code/03-pol.R", fil, foo
    output:	pol
    log: 	"logs/pol-{sid}.Rout"
    shell: '''R CMD BATCH\\
    --no-restore --no-save "--args wcs={wildcards}\
    {input[1]} {input[2]} {output}" {input[0]} {log}'''

# regions
def foo(wildcards): return(expand(
	"imgs/{did}/shapes/{sid}/{roi}",
	did=MD["did"][MD["sid"] == wildcards.sid],
	sid=wildcards.sid, roi=ROI[wildcards.sid]))
rule roi:
	input:	"code/03-roi.R", fil, x = foo
	params: lambda wc, input: pool(input.x)
	output:	roi
	log:	"logs/roi-{sid}.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {params} {output[0]}" {input[0]} {log}'''

# signatures
rule sig:
    priority: 97
    threads: 10
    input:  "code/03-sig.R", fil, 
            x = expand("meta/sig/sig-{sub}.json", sub=SUB)
    params:	lambda wc, input: pool(input.x)
    output:	sig
    log:    "logs/sig-{sid}.Rout"
    shell: '''R CMD BATCH\\
    --no-restore --no-save "--args wcs={wildcards}\
    {input[1]} {params} {output} ths={threads}" {input[0]} {log}'''


# communication
rule ccc:
	threads: 33
	input:	"code/04-ccc.R", fil
	output:	ccc
	log:    "logs/ccc-{sid}.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {output[0]} ths={threads}" {input[0]} {log}'''

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

# metastasis
rule met:
    priority: 97
    input:  "code/03-met.R",
            "outs/fil-241.rds"
    output: "outs/lv1-241.rds"
    log:    "logs/met.Rout"
    shell: '''R CMD BATCH\\
    --no-restore --no-save "--args wcs={wildcards}\
    {input[1]} {output}" {input[0]} {log}'''
rule net:
    priority: 97
    input:  "code/04-net.R",
            "outs/fil-241.rds",
            "outs/lv1-241.rds",
            "outs/lv2-242,epi.rds"
    output: "outs/lv2-241,epi.rds"
    log:    "logs/net.Rout"
    shell: '''R CMD BATCH\\
    --no-restore --no-save "--args wcs={wildcards}\
    {input[1]} {input[2]} {input[3]}\
    {output}" {input[0]} {log}'''

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
    input:	"code/06-kst.R", 
            rules.roi.output,
            rules.lv1.output,
            rules.ref.output, 
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
	params:	lambda wc, input: pool(input.x),
			lambda wc, input: pool(input.y)
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
	params:	lambda wc, input: pool(input.x),
			lambda wc, input: pool(input.y)
	output:	qbs
	log:    "logs/qbs-{sub}.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{params} {output} ths={threads}" {input[0]} {log}'''

one = expand("outs/sub-{sid},epi.rds", sid=SID)
two = expand("outs/kst-{sid}.rds", sid=SID)
rule qbt:
	threads: 30
	priority: 95
	input:	"code/00-pbs.R", x=one, y=two
	params:	lambda wc, input: pool(input.x),
			lambda wc, input: pool(input.y)
	output:	qbt
	log:    "logs/qbt.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{params} {output} ths={threads}" {input[0]} {log}'''

one = expand("outs/sub-{sid},{{sub}}.rds", sid=SID)
two = expand("outs/lv2-{sid},{{sub}}.rds", sid=SID)
rule mgs:
	threads: 30
	priority: 95
	input:	"code/00-mgs.R", x=one, y=two
	params:	lambda wc, input: pool(input.x),
			lambda wc, input: pool(input.y)
	output:	mgs
	log:    "logs/mgs-{sub}.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{params} {output} ths={threads}" {input[0]} {log}'''

# contexts
for foo in ["ctx", "cty"]:
    rule:
        priority: 95
        name:   foo
        input:	"code/06-%s.R" % foo,
                x = expand(fil, sid=SID),
                y = expand(jst, sid=SID, sub=SUB)
        params:	lambda wc, input: pool(input.x),
                lambda wc, input: pool(input.y)
        output:	expand({"ctx": ctx, "cty": cty}[foo], sid=SID)
        log:    "logs/%s.Rout" % foo
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
	input:	"code/06-trj.R", rep, roi
	output:	trj
	log:    "logs/trj-{sid}.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {input[2]} {output}" {input[0]} {log}'''	

# reprocessing
rule req:
	threads: 20
	priority: 94
	input:	"code/07-req.R", 
            rules.fil.output,
            rules.kst.output
	output:	req
	log:    "logs/req-{sid}.Rout"
	shell: '''R CMD BATCH\\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {input[2]} {output} ths={threads}" {input[0]} {log}'''

# collection
rule res:
    priority: 90
    input:  "code/09-res.R", "outs/raw-{sid}",
            expand("outs/{out}-{{sid}}.rds", 
                out=["fil", "roi", "sig", "ccc", "cty", 
                    "trj", "pro", "rep", "ist", "lv1"]),
            expand("outs/{out}-{{sid}},{sub}.rds", 
                out=["jst", "lv2"], sub=SUB)
    output: directory("outs/res-{sid}")
    log:    "logs/res-{sid}.Rout"
    shell: '''R CMD BATCH\\
    --no-restore --no-save "--args wcs={wildcards}\
    {input} {output}" {input[0]} {log}'''

# plotting =========================================

rule plt_one:
	priority: 99
	input:	"10-plt__one-{out},{plt}.R", "outs/{out}.rds"
	log:	"logs/plt__one-{out},{plt}.Rout"
	output: "plts/{out},{plt}.pdf"
	shell: '''R CMD BATCH\
	--no-restore --no-save "--args wcs={wildcards}\
	{input[1]} {output}" {input[0]} {log}'''

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

def out(by, n=None):
    o = "out" if n is None else "out"+str(n)
    d = {
        "sip": "outs/{"+o+"}-{x}.rds",
        "sid": "outs/{"+o+"}-{x}.rds",
        "sub": "outs/{"+o+"}-{x}.rds",
        "sid_sub": "outs/{"+o+"}-{x}.rds",
        "all_sip": expand("outs/{{"+o+"}}-{x}.rds", x=SIP),
        "all_sid": expand("outs/{{"+o+"}}-{x}.rds", x=SID),
        "all_sub": expand("outs/{{"+o+"}}-{x}.rds", x=SUB),
        "all_raw": expand("outs/{{"+o+"}}-{x}/se.rds", x=SID),
        "all_sid_one_sub": expand("outs/{{"+o+"}}-{y},{{x}}.rds", y=SID),
        "one_sid_all_sub": expand("outs/{{"+o+"}}-{{x}},{y}.rds", y=SUB),
        "all_sid_all_sub": expand("outs/{{"+o+"}}-{x},{y}.rds", x=SID, y=SUB)
    }
    return(d[by])

pat = "plt__{by}-{{out}},{{plt}}"
for by in ["sip", "sid", "sub", "all_sid_one_sub"]:
    rule:
        name:   "plt__%s" % by
        input:  expand("code/10-"+pat+".R", by=by), a = out(by)
        params: lambda wc, input: pool(input.a)
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
        params: lambda wc, input: pool(input.a)
        log:	expand("logs/"+pat+".Rout", by=by)
        output:	"plts/{out},{plt}.pdf"
        shell: '''R CMD BATCH\
        --no-restore --no-save "--args wcs={wildcards}\
        {params} {output}" {input[0]} {log}'''

pat = "plt__{by1}__{by2}-{{out1}},{{out2}},{{plt}}"
for by in [
    ["all_raw","all_sid"],
	["all_sid","all_sid"],
    ["all_sip","all_sip"],
	["all_sid","all_sid_all_sub"]]:
    rule:
        name:   "plt__%s__%s" % (by[0], by[1])
        input:  expand("code/10-"+pat+".R", by1=by[0], by2=by[1]),
                a = out(by[0], 1), b = out(by[1], 2)
        params: lambda wc, input: pool(input.a),
                lambda wc, input: pool(input.b)
        log:    expand("logs/"+pat+".Rout", by1=by[0], by2=by[1])
        output:	"plts/{out1},{out2},{plt}.pdf"
        shell: '''R CMD BATCH\
        --no-restore --no-save "--args wcs={wildcards}\
        {params} {output}" {input[0]} {log}'''

# plt by sid
pat = "plt__{by1}__{by2}-{{out1}},{{out2}},{{plt}}"
for by in [
    ["sid","sid"],
    ["sip","sip"],
    ["sid","one_sid_all_sub"],
    ["one_sid_all_sub","one_sid_all_sub"]]:
    rule:
        name:   "plt__%s__%s" % (by[0], by[1])
        input:  expand("code/10-"+pat+".R", by1=by[0], by2=by[1]),
                a = out(by[0], 1), b = out(by[1], 2)
        params: lambda wc, input: pool(input.a),
                lambda wc, input: pool(input.b)
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
        params: lambda wc, input: pool(input.a),
                lambda wc, input: pool(input.b)
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

# epithelial trajectory
def xxx(x):
	if x in ["sub", "jst", "lv2"]:
		rds = "outs/{x}-{sid},{sub}.rds"
	else:
		rds = "outs/{x}-{sid}.rds"
	return(expand(rds, x=x, sid=SID, sub=SUB))

qlt_trj_r = "code/10-qlt_trj-{x},{p}.R"
qlt_trj_p = "plts/trj,{x},{p}.pdf"

foo = glob_wildcards(qlt_trj_r)
for x,p in zip(foo.x, foo.p):
    rule:
        priority: 7
        input:  expand(qlt_trj_r, x=x, p=p), 
                a = xxx("trj"), 
                b = xxx("fil"), 
                c = xxx(x)
        log:    "logs/qlt_trj-"+x+","+p+".Rout"
        params:	lambda wc, input: pool(input.a),
                lambda wc, input: pool(input.b),
                lambda wc, input: pool(input.c)
        output:	expand(qlt_trj_p, x=x, p=p)
        shell: '''R CMD BATCH\
        --no-restore --no-save "--args wcs={wildcards}\
        {params} {output}" {input[0]} {log}'''

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
        params:	lambda wc, input: pool(input.a),
                lambda wc, input: pool(input.b),
                lambda wc, input: pool(input.c),
                lambda wc, input: pool(input.d)
        output:	expand(qlt_tcs_p, x=x, p=p)
        shell: '''R CMD BATCH\
        --no-restore --no-save "--args\
        {params} {output}" {input[0]} {log}'''
