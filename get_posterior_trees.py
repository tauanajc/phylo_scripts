#!/usr/bin/env python

'''
Created by Bruno de Medeiros and Tauana Cunha, Harvard University, 2017
https://github.com/brunoasm
https://github.com/tauanajc
Takes two or more .t files and combines them in one file, excluding burnin. Can also use for just one tree file to reduce number of trees.
If fossils present that need to be deleted from tree (for example from a Fossilized Birth-Death analysis),
the script will delete them if -f/-fossil flag is called. Must give a file with names of terminals.
Finally, .t files from MrBayes come with branch length plus clockrate, so values need to be converted to time. Use flag -t/--time for this.

python get_posterior_trees.py -b 25 -t -f fossil_names.txt mrbayes.run1.t mrbayes.run2.t
'''

import argparse, dendropy, math, gzip, warnings, subprocess, sys

def brlen2time(tree):
	clockrate = float(tree.annotations.get_value('clockrate'))
	for node in tree.postorder_node_iter():
		if node.edge_length:
			node.edge_length = node.edge_length / clockrate
	return tree



#to make program run faster and save memory, only trees to be kept are read into memory, the rest are not parsed
#we also assume that all files have the same number of trees
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='This program takes the output of MrBayes trees produced with FBD, removes fossils and burn-in and thins the result to a desired number of trees. Retained trees are saved in a file named posterior.trees. Optionally, TreeAnnotator or sumtrees.py can be run to obtain a summary tree from this posterior.')
    parser.add_argument('input',help = 'path to tree files (can be gzipped)', nargs = '+')
    parser.add_argument('-b', '--burnin', default = 50, type = int, help = 'Percentage of trees to discard as burnin')
    parser.add_argument('-k', '--keep', default = 10000, type = int, help = 'Maximum number of post-burnin trees to keep')
    parser.add_argument('-f','--fossil', type = str, help = 'File with fossil names to remove')
    parser.add_argument('-t','--time', default = False, action = 'store_true', help = 'Use this flag if you want to convert branch lengths to time (usually needed if trees produced in MrBayes)')
    parser.add_argument('-s','--sumtrees', default = False, action = 'store_true', help = 'summarizes posterior trees in an MCC tree using sumtrees.py')
    parser.add_argument('-a','--treeannotator', type = str, help = 'path to beast.jar if summarization with treeannotator is desired')
    
    args = parser.parse_args()
    #args =  parser.parse_args('-s -t --fossil fossils.txt -k 10 full.nex.run1.t'.split())
    
    #read fossils file if present
    if args.fossil:
        fossil_labels = [line.rstrip('\n') for line in open(args.fossil, 'r')]
    
    print 'CHOOSING TREES TO BE READ'
    sys.stdout.flush()
    
    #first, read only first tree to make a taxon name space
    ntrees_per_file = 0 #records number of trees
    to_read = True
    
    lines = []
    if args.input[0][-3:] == '.gz':
        infile = gzip.open(args.input[0], 'rb')
    else:
        infile = open(args.input[0], 'r')
    
    for line in infile:
        if to_read: lines.append(line) #keeps first block of nexus plus first tree only
        if 'tree gen' in line or 'tree STATE' in line:
            ntrees_per_file += 1
            to_read = False
    infile.close()
    
    first = dendropy.Tree.get(data = '\n'.join(lines), schema = 'nexus', preserve_underscores=True) #gets first tree
    taxa = first.taxon_namespace #gets taxa names from first tree
    
    
    #now, calculate which trees to keep
    nfiles = len(args.input) #number of input files
    n_keep_per_file = int(math.ceil(args.keep/float(nfiles))) #k/(number of input files)
    
    burnin = int(math.ceil(args.burnin/100.0 * ntrees_per_file)) #number of trees to eliminate as burnin
    thinning = max(1,int(math.floor((float(ntrees_per_file) - burnin) / n_keep_per_file))) #to get sample tree every thinning# trees
    #if thinning==0: thinning = thinning + 1 #for when math.floor leads to thinning = 0
    
    #if number of trees to keep is multiple of number of files, this results in the exact number of trees, with at least the specified burnin
    #if not, will result in a few more trees than requested
    lines_to_keep = [tree_idx + len(lines) - 1 for tree_idx in range(burnin,ntrees_per_file,thinning)[-n_keep_per_file:]]
    				#tree_idx is iterator on tree line number
 							#(-1) to remove first tree in nexus block
 											#range = seq starting in burnin, going to total number of trees, skipping thinning# trees
    
    
    #read trees and save to a dendropy tree list
    print 'READING TREES AND REMOVING FOSSILS'
    sys.stdout.flush()
    nread = 0
    
    treelist = dendropy.TreeList()
    for infilepath in args.input:
        if infilepath[-3:] == '.gz':
            infile = gzip.open(infilepath, 'rb')
        else:
            infile = open(infilepath, 'r')
    
   	#enumerate iterates in a file, returning a tuple with two values: 1) the index of the line and 2) the string in the line
        for k,line in enumerate(infile):
            if k in lines_to_keep: #when reads line number that matches lines_to_keep, appends tree to output file after the nexus block
                nread += 1
                sys.stdout.write('TREE ' + '{:2d}'.format(nread) + '\r')
                sys.stdout.flush()
    
                tree = dendropy.Tree.get(data = '\n'.join(lines[:-1] + [line]), schema = 'nexus', taxon_namespace = taxa, preserve_underscores=True)
                
                if args.time:
                    tree = brlen2time(tree)
                    #treelist.append(tree)            
                
                if args.fossil: #remove fossils
                    #since we will update bipartitions and suppress_unifurcations, we need to add a fake edge to the seed node to avoid it disappearing
                    tree.seed_node.add_child(dendropy.Node(taxon=dendropy.Taxon('root'),edge_length=1))
                    tree.update_taxon_namespace()
                    #now we remove fossils
                    #print fossil_labels
                    tree.prune_taxa_with_labels(fossil_labels, update_bipartitions=True, suppress_unifurcations=True)
                    #and remove the fake edge
                    tree.prune_taxa_with_labels(['root'],update_bipartitions=True,suppress_unifurcations=False)
                    #finally, remove all non-used taxa from the taxon namespace
                    tree.purge_taxon_namespace()
    
                    
                    #sometimes this causes rounding errors and downstream applications consider the tree non-ultrametric
                    #first, test if ultrametric with precision 0.001
                    try:
                        tree.calc_node_ages(ultrametricity_precision=1e-03)
                    except ValueError as err:
                        print err
                        warnings.warn('Some tree is not ultrametric!')
                    else:    
                        #if it is, let's prune a little bit from each tip to make it ultrametric
                        tip_dists = {leaf:leaf.distance_from_root() for leaf in tree.leaf_node_iter()}
                        min_tip_dist = min([i for k, i in tip_dists.iteritems()])
                        for node, this_tip_dist in tip_dists.iteritems():
                            node.edge.length = node.edge.length - (this_tip_dist - min_tip_dist)
                    #tree.print_plot()
                            
                treelist.append(tree)
                            
                    
    infile.close()
    
    treelist.purge_taxon_namespace()
    
    #save nexus
    print 'SAVING CHOSEN TREES'
    sys.stdout.flush()
    treelist.write(path = 'posterior.trees', schema = 'nexus', 
                translate_tree_taxa = True, 
                unquoted_underscores = True, 
                ignore_unrecognized_keyword_arguments = True,
                suppress_taxa_block = True
                )
    
    del treelist
    
    #summarize
    if args.sumtrees:
        print 'GETTING MCC TREE WITH SUMTREES.PY'
        run_sumtrees = subprocess.Popen('sumtrees.py --summary-target=mcct --edges=median-age -o sumtrees_mcc.tre posterior.trees'.split())
        run_sumtrees.communicate()
        
        print 'ADJUSTING ANNOTATIONS FOR FIGTREE'
        sys.stdout.flush()
        tree = dendropy.Tree.get(path = 'sumtrees_mcc.tre', schema = 'nexus')
        tree.write(path = 'sumtrees_mcc.tre', schema = 'nexus', 
                translate_tree_taxa = True, 
                unquoted_underscores = True, 
                suppress_taxa_block = False)
                
    if args.treeannotator:
        print 'GETTING MCC TREE WITH TREEANNOTATOR'
        treean_command = 'beast.app.treeannotator.TreeAnnotator -burnin 0 -heights median posterior.trees treeannotator_mcc.tre'.split()
        beast_command = ['java','-cp',args.treeannotator]
        beast_command.extend(treean_command)
        run_sumtrees = subprocess.Popen(beast_command)
        run_sumtrees.communicate()

        
    print 'DONE'
