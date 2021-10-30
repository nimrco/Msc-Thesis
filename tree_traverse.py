import subprocess


'''
tree = read('tree_file', format='newick', into=TreeNode)
for node in tree.postorder():
    print(node.name)
'''

'''
tree_test = TreeNode.read(["((a,b,(c,d)e)f,(g,h)i)root;"])
print(tree_test.ascii_art())
f = lambda n: [n.name] if n.is_tip() else []
tree_test.cache_attr(f, 'tip_names')
'''

res = subprocess.run("Rscript pack.R", shell=True)
print("stdout: {}".format(res.stdout))
print("stderr: {}".format(res.stderr))
