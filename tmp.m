     % Create a figure window
      figure, hold on;
     % x,y,z line coordinates
      x=60*sind(0:20:360); y=60*cosd(0:20:360); z=60*cosd((0:20:360)*2); 
     % plot the first line
      h=plot3t(x,y,z,linspace(1,5,length(x)));
     % Shade the line
      set(h, 'FaceLighting','phong','SpecularColorReflectance', 0, 'SpecularExponent', 50, 'DiffuseStrength', 1);
     % Set figure rotation
      view(3); axis equal;
     % Set the material to shiny and enable light
      material shiny; camlight;
     