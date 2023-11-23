function scalebar(xc,yc,zc, w, label)
  x = [xc, xc-w, nan, xc, xc  , nan, xc, xc  ];
  y = [yc, yc  , nan, yc, yc-w, nan, yc, yc  ];
  z = [zc, zc  , nan, zc, zc  , nan, zc, zc+w];    
  hl = line(x,y,z);
  hl.LineWidth = 2;
  ht = text(xc,yc,zc,[num2str(w), ' ', label]);
  ht.FontSize = 18;
  ht.Color = hl.Color;
  ht.VerticalAlignment = 'bottom';
  