for file in /home/anacarolina/MLP/instances/*
        do
          ./tsp "$file" >> testesMLP.txt
        done 