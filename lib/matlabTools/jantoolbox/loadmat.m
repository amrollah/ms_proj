function A = loadmat(file)
AA = load(file);
A = sparse(AA(2:end,1)+1,AA(2:end,2)+1,AA(2:end,3),AA(1,1),AA(1,2));
