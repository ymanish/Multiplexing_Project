import os
import re
path = "."
unsucessful_file = open('run_again_org.txt', 'w')

for root,d_names,f_names in os.walk(path):
    print (f_names)
    count=1
    for i in f_names:
        # print(i)
        if re.search("out$", i):
            f = open(i, 'r')
            last_line = f.readlines()[-1]
            if re.search("^Combine_Charts", last_line):
                print('ok')
                f.close()

            else:

                f = open(i, 'r')
                tobe_down = f.readlines()[4].split(' ')[2]
                print(tobe_down)
                unsucessful_file.write(str(count))
                unsucessful_file.write(',')
                unsucessful_file.write(tobe_down)
                unsucessful_file.write('\n')
                f.close()
                count=count+1

