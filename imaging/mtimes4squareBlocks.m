function out = mtimes4squareBlocks( O, w, mask )

x = w(1:end/2);
y = w(1+end/2:end);
out = [...
    maskedMtimes(O{1,1},x,mask) + maskedMtimes(O{1,2},y,mask)
    maskedMtimes(O{2,1},x,mask) + maskedMtimes(O{2,2},y,mask)
    ];

end

