function out = mtimes4squareBlocks( O, w )

x = w(1:end/2);
y = w(1+end/2:end);
out = [...
    O{1,1}*x + O{1,2}*y
    O{2,1}*x + O{2,2}*y
    ];

end

