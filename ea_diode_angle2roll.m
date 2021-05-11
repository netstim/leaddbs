%% converts orientation angle from CT slice to roll angle of the lead according to yaw and pitch
function roll = ea_diode_angle2roll(angle,yaw,pitch)
        roll = (sin(angle) * cos(pitch)) / ((cos(angle) * cos(yaw)) - (sin(angle) * sin(yaw) * sin(pitch)));  % see Sitz et al. 2017
        roll = atan(roll);        
        if angle < pi && roll < 0 && angle - roll > pi/2
            roll = roll + pi;
        end
        if angle > pi && roll > 0 && angle - roll > pi/2
            roll = roll - pi;
        end
end