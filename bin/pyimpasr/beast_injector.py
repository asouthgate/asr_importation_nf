from Bio.SeqIO.FastaIO import SimpleFastaParser
import sys
import pandas as pd
import datetime as dt


def get_clock_rate(clock_rate):
    return("""<parameter id="default.clock.rate" value="{cr}" lower="0.0"/>""".format(cr=clock_rate))


def get_taxon_block(name,date,loc,part):
    return ("""        <taxon id="{n}">
            <date value="{d}" direction="forwards" units="years"/>
            <attr name="{p}_loc">
                {l}
            </attr>
        </taxon>""".format(n=name,d=date,l=loc,p=part))

get_taxon_line = lambda name: """        <taxon idref="{n}"/>""".format(n=name)

def get_taxon_dates_and_locs(taxa_data):
    lines = []
    for taxon, date, loc, part in taxa_data:
        lines.append(get_taxon_block(taxon, date, loc,part))
    return "\n".join(lines)        

def get_partition(partition_name, taxon_names):
    lines = ["""    <taxa id="{pn}.taxa">""".format(pn=partition_name)]
    for tn in taxon_names: lines.append(get_taxon_line(tn))
    lines.append("""    </taxa>""")
    return "\n".join(lines)

def get_first_tree_model(partition_names):
    return """        <treeModel idref="{p}.treeModel"/>""".format(p=partition_names[0])

def get_partitions(partition_names, taxon_sets):
    lines = []
    for pi, pn in enumerate(partition_names):
        lines.append(get_partition(pn, taxon_sets[pi]))
    return "\n".join(lines)        

def get_aln_name(partition_name):
    return partition_name + ".aln"

def get_alignment(partition_name, taxon_records):
    aln_name = get_aln_name(partition_name)
    lines = ["""    <alignment id="{an}" dataType="nucleotide">""".format(an=aln_name)]
    for header, seq in taxon_records:
        lines.append("""        <sequence>""")
        lines.append("""            <taxon idref="{h}"/>""".format(h=header))
        lines.append("""            {s}""".format(s=seq))
        lines.append("""        </sequence>""")
    lines.append("""    </alignment>""")
    return "\n".join(lines)

def get_alignments(partition_names, record_lists):
    return "\n".join(get_alignment(pn, record_lists[pi]) for pi,pn in enumerate(partition_names))

def get_single_alignment(records):
    return get_alignment("alignment", records)


def get_patterns(partition_names):
    lines = []
    for pn in partition_names:
        patt_name = pn + ".pattern"
        aln_name = get_aln_name(pn)
        lines.append("""    <patterns id="{pattern_name}" from="1" strip="false">
        <alignment idref="{an}"/>
    </patterns>\n""".format(pattern_name=patt_name, an=aln_name))
    return "\n".join(lines)

def get_attrib_pattern_name(pn): return pn + "_loc.pattern"

def get_attrib_pattern_dt(pn): return "first.dataType"

def get_taxa_set(pn): return pn + ".taxa"

def get_starting_tree_name(pn): return pn + ".startingTree"

def get_attribute_patterns(partition_names):
    lines = []
    for pn in partition_names:
        attrib_pattern_name = get_attrib_pattern_name(pn)
        attrib_pattern_dt = get_attrib_pattern_dt(pn)
        taxa_set = get_taxa_set(pn)
        lines.append("""    <attributePatterns id="{apn}" attribute="{pn}_loc">""".format(apn=attrib_pattern_name, pn=pn))
        lines.append("""        <taxa idref="{taxnames}"/>""".format(taxnames=taxa_set))
        lines.append("""        <generalDataType idref="{apdt}"/>""".format(apdt=attrib_pattern_dt))
        lines.append("""    </attributePatterns>\n""")
    return "\n".join(lines)

def get_starting_trees(partition_names):
    lines = []
    for pn in partition_names:
        starting_tree_name = get_starting_tree_name(pn)
        stn = pn + ".startingTree"
        tax = pn + ".taxa"
        lines.append("""    <coalescentSimulator id="{stn}" height="0.9">""".format(stn=starting_tree_name))
        lines.append("""        <taxa idref="{tax}"/>""".format(tax=tax))
        lines.append("""        <constantSize idref="initialDemo"/>""")
        lines.append("""    </coalescentSimulator>\n""")
    return "\n".join(lines)

def get_tree_updown(partition_names):
    lines = []
    for pn in partition_names:
        lines.append("""			<up>
				<parameter idref="{p}.treeModel.allInternalNodeHeights"/>
			</up>""".format(p=pn))
    return "\n".join(lines)

def get_first_treemodel_reference(partition_names):
    return """                    <treeModel idref="{p}.treeModel"/>""".format(p=partition_names[0])

def get_tree_models_and_statistics(partition_names):
    lines = []
    for pn in partition_names:
        lines.append("""    <treeModel id="{p}.treeModel">
        <coalescentTree idref="{p}.startingTree"/>
        <rootHeight>
            <parameter id="{p}.treeModel.rootHeight"/>
        </rootHeight>
        <nodeHeights internalNodes="true">
            <parameter id="{p}.treeModel.internalNodeHeights"/>
        </nodeHeights>
        <nodeHeights internalNodes="true" rootNode="true">
            <parameter id="{p}.treeModel.allInternalNodeHeights"/>
        </nodeHeights>
    </treeModel>\n""".format(p=pn))
        lines.append("""    <treeLengthStatistic id="{p}.treeLength">
        <treeModel idref="{p}.treeModel"/>
    </treeLengthStatistic>\n""".format(p=pn))
        lines.append("""    <tmrcaStatistic id="{p}.age(root)" absolute="true">
        <treeModel idref="{p}.treeModel"/>
    </tmrcaStatistic>\n""".format(p=pn))
    return "\n".join(lines)

def get_population_tree_within_gmrf_skygrid(max_pn):
    return """            <treeModel idref="{pn}.treeModel"/>""".format(pn=max_pn) 


def get_likelihood_for_tree_given_sequence_data(partition_names):
    lines = []
    for pn in partition_names:
        lines.append("""    <treeDataLikelihood id="{p}.treeLikelihood" useAmbiguities="false">
        <partition>
            <patterns idref="{p}.pattern"/>
            <siteModel idref="siteModel"/>
        </partition>
        <treeModel idref="{p}.treeModel"/>
        <strictClockBranchRates idref="UK1018.part.m.branchRates"/>
    </treeDataLikelihood>\n""".format(p=pn))
    return "\n".join(lines)



def get_ancestral_tree_likelihood(partition_names):
    lines = []
    for pn in partition_names:
        lines.append("""    <ancestralTreeLikelihood id="{p}_loc.treeLikelihood" stateTagName="{p}.states" useUniformization="true" saveCompleteHistory="false" logCompleteHistory="false">
        <attributePatterns idref="{p}_loc.pattern"/>
        <treeModel idref="{p}.treeModel"/>
        <siteModel idref="first.siteModel"/>
        <generalSubstitutionModel idref="first.model"/>
        <generalSubstitutionModel idref="first.model"/>
        <rateEpochBranchRates idref="epochClockRates"/>
    </ancestralTreeLikelihood>\n""".format(p=pn))
    return "\n".join(lines)

def get_tree_operators(partition_names):
    lines = []
    for pn in partition_names:
        lines.append("""        <subtreeSlide size="1.0" gaussian="true" weight="30">
            <treeModel idref="{p}.treeModel"/>
        </subtreeSlide>
        <narrowExchange weight="30">
            <treeModel idref="{p}.treeModel"/>
        </narrowExchange>
        <wideExchange weight="3">
            <treeModel idref="{p}.treeModel"/>
        </wideExchange>
        <wilsonBalding weight="3">
            <treeModel idref="{p}.treeModel"/>
        </wilsonBalding>
        <scaleOperator scaleFactor="0.75" weight="3">
            <parameter idref="{p}.treeModel.rootHeight"/>
        </scaleOperator>
        <uniformOperator weight="30">
            <parameter idref="{p}.treeModel.internalNodeHeights"/>
        </uniformOperator>\n""".format(p=pn))
    return "\n".join(lines)

def get_tree_root_height_priors(partition_names):
    lines = []
    for pn in partition_names:
        lines.append("""                <uniformPrior lower="0.0" upper="1.0">
                    <parameter idref="{p}.treeModel.rootHeight"/>
                </uniformPrior>\n""".format(p=pn))
    return "\n".join(lines)
    
def get_tree_likelihood_references(partition_names):
    lines = []
    for pn in partition_names:
        lines.append("""                <treeDataLikelihood idref="{p}.treeLikelihood"/>""".format(p=pn))
    return "\n".join(lines)

def get_ancestral_tree_likelihood_references(partition_names):
    lines = []
    for pn in partition_names:
        lines.append("""                <ancestralTreeLikelihood idref="{p}_loc.treeLikelihood"/>""".format(p=pn))
    return "\n".join(lines)

def get_tree_rootheight_logs(partition_names):
    lines = []
    for pn in partition_names:
        lines.append("""            <parameter idref="{p}.treeModel.rootHeight"/>
            <tmrcaStatistic idref="{p}.age(root)"/>
            <treeLengthStatistic idref="{p}.treeLength"/>\n""".format(p=pn))
    return "\n".join(lines)

def get_age_root_logs(partition_names):
    lines = []
    for pn in partition_names:
        lines.append("""            <column label="{p}.age(root)" sf="6" width="12">
                <tmrcaStatistic idref="{p}.age(root)"/>
            </column>""".format(p=pn))
    return "\n".join(lines)

def get_tree_data_likelihood_log_ref(partition_names):
    lines = []
    for pn in partition_names:
        lines.append("""            <treeDataLikelihood idref="{p}.treeLikelihood"/>""".format(p=pn))
    return "\n".join(lines)

def get_ancestral_likelihood_log_ref(partition_names):
    lines = []
    for pn in partition_names:
        lines.append("""            <ancestralTreeLikelihood idref="{p}_loc.treeLikelihood"/>\n""".format(p=pn))
    return "\n".join(lines)

def get_fixed_tree(tree):
    s = """
                <newick id="startingTree" usingDates="false" units="time" usingHeights="true">
                    {nwk}
                </newick>""".format(nwk=tree)
    return s


def get_tree_logs(partition_names):
    lines = []
    for pn in partition_names:
        lines.append("""        <logTree id="{p}.treeFileLog" logEvery="1000" nexusFormat="true" fileName="{p}.trees" sortTranslationTable="true">
            <treeModel idref="{p}.treeModel"/>
            <trait name="rate" tag="UK1018.part.m.rate">
                <strictClockBranchRates idref="UK1018.part.m.branchRates"/>
            </trait>
<!--
            <trait name="rate" tag="first.rate">
                <strictClockBranchRates idref="first.branchRates"/>
            </trait>
-->
            <joint idref="joint"/>
            <trait name="{p}.states" tag="first">
                <ancestralTreeLikelihood idref="{p}_loc.treeLikelihood"/>
            </trait>
        </logTree>\n""".format(p=pn))
    return "\n".join(lines)

funcd = {"!@TAXON_DATES_AND_LOCS@!":get_taxon_dates_and_locs,
        "!@PARTITIONS@!":get_partitions,
        "!@FIRST_TREEMODEL_REFERENCE@!":get_first_treemodel_reference,
        "!@ALIGNMENTS@!":get_alignments,
        "!@SINGLE_ALIGNMENT@!":get_single_alignment,
        "!@PATTERNS@!":get_patterns,
        "!@ATTRIBUTE_PATTERNS@!":get_attribute_patterns,
        "!@STARTING_TREES@!":get_starting_trees,
        "!@TREE_MODELS_AND_STATISTICS@!":get_tree_models_and_statistics,
        "!@POPULATION_TREE_WITHIN_GMRF_SKYGRID@!":get_population_tree_within_gmrf_skygrid,
        "!@LIKELIHOOD_FOR_TREE_GIVEN_SEQUENCE_DATA@!":get_likelihood_for_tree_given_sequence_data,
        "!@ANCESTRAL_TREE_LIKELIHOOD@!":get_ancestral_tree_likelihood,
        "!@TREE_OPERATORS@!":get_tree_operators,
        "!@TREE_ROOT_HEIGHT_PRIORS@!":get_tree_root_height_priors,
        "!@TREE_LIKELIHOOD_REFERENCES@!":get_tree_likelihood_references,
        "!@ANCESTRAL_TREE_LIKELIHOOD_REFERENCES@!":get_ancestral_tree_likelihood_references,
        "!@TREE_LOGS@!":get_tree_logs,
        "!@TREE_DATA_LIKELIHOOD_LOG_REF@!":get_tree_data_likelihood_log_ref,
        "!@ANCESTRAL_LIKELIHOOD_LOG_REF@!":get_ancestral_likelihood_log_ref,
        "!@TREE_ROOTHEIGHT_LOGS@!":get_tree_rootheight_logs,
        "!@TREE_LOGS@!":get_tree_logs,
        "!@FIRST_TREE_MODEL@!":get_first_tree_model,
        "!@TREE_UPDOWN@!":get_tree_updown,
        "!@AGE_ROOT_LOGS@!":get_age_root_logs,
        "!@FIXED_TREE@!":get_fixed_tree,
        "!@CLOCK_RATE@!":get_clock_rate

}

if __name__ == "__main__":
    import argparse as ap
    from Bio import Phylo
    args.add_argument("--fixed_tree", action="store_true", default=False)
    args = ap.parse_args()
    if args.fixed_tree:
        assert ftreefile.endswith("newick")
        with open(args.fixed_tree) as f:
            fixed_tree = [l.rstrip("\n") for l in f][0]
            
    partition_folder_names = [fn for fn in sys.argv[1:]]
    partition_names = [fn.split("/")[-1] for fn in partition_folder_names]
    print(partition_folder_names)
    for pfn in partition_folder_names: assert "." not in pfn
    record_lists = []
    taxon_sets = []
    for pfn in partition_folder_names:
        pref = pfn.split("/")[-1]
        fn = pfn + "/" + pref + ".m.fa"
        assert fn.endswith(".fa")
        with open(fn) as f:
            recs = [tup for tup in SimpleFastaParser(f)]
            record_lists.append(recs)
            taxon_sets.append([h for h,s in recs])
    taxa_data = []
    start = dt.datetime.strptime("2020-01-01", "%Y-%m-%d") 
    for pfn in partition_folder_names:
        pref = pfn.split("/")[-1]
        csvfn = pfn + "/" + pref + ".csv"
        with open(csvfn) as f:
            lines = [l.strip() for l in f][1:]
            for l in lines:
                tax,date,loc,part = l.strip().split(",")
                date = 2020 + (dt.datetime.strptime(date, "%Y-%m-%d")-start).days/365.0
                taxa_data.append((tax,date,loc,part))
    
    max_partition_name = max([z for z in zip(partition_names, taxon_sets)])[0]
    with open("scripts/template_epoch.xml") as f:
        for l in f:
            l = l.rstrip()
            if l.startswith("!@"):
                func = funcd[l.strip()]
                if l.startswith("!@TAXON_DATES_AND_LOCS@!"):
                    print(func(taxa_data))
                elif l.startswith("!@PARTITIONS@!"):
                    print(func(partition_names, taxon_sets))
                elif l.startswith("!@ALIGNMENTS@!"):
                    print(func(partition_names, record_lists))
                elif l.startswith("!@POPULATION_TREE_WITHIN_GMRF_SKYGRID@!"):
                    print(func(max_partition_name))
                elif l.startswith("!@FIXED_TREE@!"):
                    print(func(fixed_tree))
                else:
                    print(func(partition_names))
            else:
                print(l)
