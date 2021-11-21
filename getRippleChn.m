%% Convert from Utah channel to Ripple channel
function chn = getRippleChn(c)

% Get bank and channel within bank 
rem = mod((c-1),32)+1;
bank = (c-rem)/32;

% Each connector on Ripple can have up to 128 channels
chn = bank*128 + rem;

end