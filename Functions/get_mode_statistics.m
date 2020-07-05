


function [out, c_mag, c_phase] = get_mode_statistics(c, dim)

%{
    This function calculates the mean and standard deviation in mode
    magnitudes (energy) and phase. The code uses the input "dim" to
    determine what dimension statistics are being generated from. 


    Note: Averaging the coefficients then taking the magnitude is not the
    same as taking the magnitude and then averaging the coefficients. Both 
    the magnitude and phase operators are nonlinear. The correct method is 
    first calculate the magnitude and phase then generate statistics. This
    prevents coherent effects.

%}
c_mag = abs(c).^2;
c_mean_mag = mean(c_mag, dim);
c_std_mag = std(c_mag, [], dim);

c_phase = (180/pi).*angle(c);
c_mean_phase = mean(c_phase, dim);
c_std_phase = std(c_phase, [], dim);

out.mean.mag = (c_mean_mag);
out.mean.phase = (c_mean_phase);
out.std.mag = (c_std_mag);
out.std.phase = (c_std_phase);

end


























