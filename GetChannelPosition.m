function [column, row, electrode] = GetChannelPosition(channel)
% This function takes the electrode channel and returns the row and column for a
% subplot of 10x10 plots, along with the electrode number
%
%//   legend
% //   col - 0 based column from left to right
% //   row - 0 based row from bottom to top
% //   bank - bank name - values can be A=1 B=2 C=3 or D
% //   elec - 1 based electrode number within the bank - values can be 1-32
% //   label - label used to rename channels in Central (optional)
%
% andrew 10 may 2016
%
%       C   R  Bank E   Ch
map = [ 1	2	3	1	96;
        1	3	3	3	95;
        1	4	3	5	94;
        1	5	3	7	93;
        1	6	3	9	92;
        1	7	3	11	91;
        1	8	3	13	90;
        1	9	3	15	89;
        2	1	1	2	88;
        2	2	3	2	87;
        2	3	3	4	86;
        2	4	3	6	85;
        2	5	3	8	84;
        2	6	3	10	83;
        2	7	3	12	82;
        2	8	3	14	81;
        2	9	3	16	80;
        2	10	3	17	79;
        3	1	1	1	78;
        3	2	2	1	77;
        3	3	2	3	76;
        3	4	2	5	75;
        3	5	2	7	74;
        3	6	2	9	73;
        3	7	2	13	72;
        3	8	3	18	71;
        3	9	3	20	70;
        3	10	3	19	69;
        4	1	1	3	68;
        4	2	2	2	67;
        4	3	2	4	66;
        4	4	2	6	65;
        4	5	2	8	64;
        4	6	2	11	63;
        4	7	2	15	62;
        4	8	2	17	61;
        4	9	3	22	60;
        4	10	3	21	59;
        5	1	1	4	58;
        5	2	1	7	57;
        5	3	1	5	56;
        5	4	2	16	55;
        5	5	2	10	54;
        5	6	2	12	53;
        5	7	2	19	52;
        5	8	2	21	51;
        5	9	3	23	50;
        5	10	3	24	49;
        6	1	1	6	48;
        6	2	1	9	47;
        6	3	1	17	46;
        6	4	1	15	45;
        6	5	2	18	44;
        6	6	2	14	43;
        6	7	2	24	42;
        6	8	2	23	41;
        6	9	3	25	40;
        6	10	3	26	39;
        7	1	1	8	38;
        7	2	1	11	37;
        7	3	1	13	36;
        7	4	1	19	35;
        7	5	2	22	34;
        7	6	2	20	33;
        7	7	2	26	32;
        7	8	2	25	31;
        7	9	3	27	30;
        7	10	3	28	29;
        8	1	1	10	28;
        8	2	1	12	27;
        8	3	1	23	26;
        8	4	1	25	25;
        8	5	1	21	24;
        8	6	2	30	23;
        8	7	2	28	22;
        8	8	2	27	21;
        8	9	3	30	20;
        8	10	3	29	19;
        9	1	1	14	18;
        9	2	1	16	17;
        9	3	1	20	16;
        9	4	1	27	15;
        9	5	1	29	14;
        9	6	1	31	13;
        9	7	2	32	12;
        9	8	2	29	11;
        9	9	2	31	10;
        9	10	3	32	9;
        10	2	1	18	8;
        10	3	1	22	7;
        10	4	1	24	6;
        10	5	1	26	5;
        10	6	1	28	4;
        10	7	1	30	3;
        10	8	1	32	2;
        10	9	3	31	1];
    
map(:,4) = map(:,4)+(map(:,3)-1)*32; % convert the channels to overall
    
inds = find(map(:,4)==channel);
column = map(inds,1);
row = map(inds,2);
electrode = map(inds,4);