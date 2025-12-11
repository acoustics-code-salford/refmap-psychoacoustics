function phon = shmSone2Phon(sone)
% phon = shmSone2Phon(sone)
%
% Returns the loudness level(s) in phon corresponding with sensory loudness
% value(s) in sone, according to the Sottek Hearing Model as defined in
% ECMA-418-2:2025.
%
% Inputs
% ------
% sone : float or array of floats
%   Loudness value(s), which must always be positive.
%
% Returns
% -------
% 
% phon : float or array of floats
%   Loudness level(s) corresponding with input sone.
%
% Assumptions
% -----------
% The input sone value(s) represents loudness according to the Sottek
% Hearing Model. The output phon is returned based on linear
% interpolation between sone-to-phon values calculated at discrete
% free-field sound pressure levels. Values outside the discrete range are
% extrapolated.
%
% Requirements
% ------------
% None
%
% Ownership and Quality Assurance
% -------------------------------
% Author: Mike JB Lotinga (m.j.lotinga@edu.salford.ac.uk)
% Institution: University of Salford
%
% Date created: 08/12/2025
% Date last modified: 11/12/2025
% MATLAB version: 2023b
%
% Copyright statement: This file and code is part of work undertaken within
% the RefMap project (www.refmap.eu), and is subject to licence as detailed
% in the code repository
% (https://github.com/acoustics-code-salford/refmap-psychoacoustics)
%
% As per the licensing information, please be aware that this code is
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%
% Checked by:
% Date last checked:
%
%% Arguments validation
    arguments (Input)
        sone (:, :) double {mustBeInRange(sone, 0, 200)}
    end

%% Define constants
phonRange = [0; 0.5; 1; 1.5; 2; 2.2; 3; 3.5; 4; 4.5; 5; 5.5; 6; 6.5; 7.5;...
             8.5; 10; 12.5; 15; 17.5; 20; 22.5; 25; 27.5; 30; 32.5; 35;...
             37.5; 40; 42.5; 45; 47.5; 50; 52.5; 55; 57.5; 60; 62.5; 65;...
             67.5; 70; 72.5; 75; 77.5; 80; 82.5; 85; 87.5; 90; 92.5; 95;...
             97.5; 100; 102.5; 105; 107.5; 110; 112.5; 115; 117.5; 120;...
             122.5; 125; 127.5; 130];

soneRange = [0.0131728008755640; 0.0147017002514385; 0.0163155244479395;...
             0.0180184921907174; 0.0198149610248851; 0.0205607307792728;...
             0.0237064956262933; 0.0259708071223978; 0.0284830012871350;...
             0.0313141063933990; 0.0342960992261228; 0.0374356436074735;...
             0.0407395216419739; 0.0442146129328367; 0.0517062900971077;...
             0.0599666456395902; 0.0739198656236416; 0.102309687878046;...
             0.137973827735358; 0.182143060528872; 0.235170666512493;...
             0.296369339360974; 0.366110735736481; 0.444665058185281;...
             0.533067637132950; 0.630859834999549; 0.739628579303990;...
             0.861665973977183; 0.999999649261453; 1.15446083738014;...
             1.32574011950862; 1.51478530574131; 1.72345292664480;...
             1.95194293008535; 2.20302763638237; 2.48230837837064;...
             2.79349557288218; 3.14141044540974; 3.53640968939006;...
             3.98960969379179; 4.50464193047483; 5.09216519681601;...
             5.76494189010932; 6.53417875936521; 7.41860357207810;...
             8.43076700943934; 9.58630438005357; 10.9139993638759;...
             12.4391485682283; 14.2012670012497; 16.2487266878483;...
             18.6386457583573; 21.4505973020449; 24.7648343591600;...
             28.6875064462480; 33.3311039332838; 38.8344811939958;...
             45.3667890880497; 53.1261734359328; 62.3737509192457;...
             73.4272890171466; 86.6550077521385; 102.511893200588;...
             121.538024632640; 144.383251212679];

if ismatrix(sone)  % flatten matrix to vector
    phonReshape = true;
    soneSize = size(sone);
    sone = sone(:);
end

phon = -Inf(size(sone));

mask = sone ~= 0;

phon(mask) = interp1(soneRange, phonRange, sone(mask), 'linear', 'extrap');

if phonReshape  % reshape matrix to vector
    phon = reshape(phon, soneSize);
end

end

