import ast

with open("exmat/ja.23.2.txt", "r") as nap_ja:
    with open("mat/ja.txt", "r") as ja:
        lines_nap = nap_ja.readlines()
        lines_ja = ja.readlines()
        
        for i in range(len(lines_nap)):
            if ast.literal_eval(lines_nap[i]) != ast.literal_eval(lines_ja[i]):
                print(i + 1)
                
                