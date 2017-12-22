function [map_base] = get_map_base(name)
%GET_MAP_BASE Returns a mapping base for use in stochastic models
%   [map_base] = GET_MAP_BASE(name)
%   returns a particular mapping base from the given name to use with
%   stochastic models.
%
%   The parameter name is a string containing the name of the mapping model
%   to be used. For now the mappings defined are:
%   - 'F9': the standard 8-bit model, each bit having a mapping. First
%   base is constant 1 to allow for arbitrary shift, second base is the most
%   significant bit, third base is the next bit, etc.
%   - 'F17': similar to F9 but for 16-bit model. First base is constant "1"
%   and then the rest are the bits, from most significant to least
%   significant, of a 16-bit value.
%   - 'F17xor': model F17 plus 8 bits corresponding to the xor between the
%   2 bytes in the 16-bit value modeled in F17.
%   - 'F17tran': model F17 plus 8 bits corresponding to transitions 0->1
%   for byte1->byte2 and transitions 1->0 in byte1->byte2. This is useful
%   only when F17 actually models two consecutive bytes rather than a full
%   16-bit value.
%
%   The map_base returned is a cell of size u x 1, where u is the number of
%   bases for the given mapping. Each base is a function of the form:
%       map_base{j} = func(v)
%   taking as parameter a value v and returning the mapping
%   of these for the base j.

%% Return base
if strcmp(name, 'F9')
    map_base = cell(9, 1);
    map_base{1} = @(v)(1);
    map_base{2} = @(v)(bitget(v, 8));
    map_base{3} = @(v)(bitget(v, 7));
    map_base{4} = @(v)(bitget(v, 6));
    map_base{5} = @(v)(bitget(v, 5));
    map_base{6} = @(v)(bitget(v, 4));
    map_base{7} = @(v)(bitget(v, 3));
    map_base{8} = @(v)(bitget(v, 2));
    map_base{9} = @(v)(bitget(v, 1));
elseif strcmp(name, 'F17')
    map_base = cell(17, 1);
    map_base{1} = @(v)(1);
    map_base{2} = @(v)(bitget(v, 16));
    map_base{3} = @(v)(bitget(v, 15));
    map_base{4} = @(v)(bitget(v, 14));
    map_base{5} = @(v)(bitget(v, 13));
    map_base{6} = @(v)(bitget(v, 12));
    map_base{7} = @(v)(bitget(v, 11));
    map_base{8} = @(v)(bitget(v, 10));
    map_base{9} = @(v)(bitget(v, 9));
    map_base{10} = @(v)(bitget(v, 8));
    map_base{11} = @(v)(bitget(v, 7));
    map_base{12} = @(v)(bitget(v, 6));
    map_base{13} = @(v)(bitget(v, 5));
    map_base{14} = @(v)(bitget(v, 4));
    map_base{15} = @(v)(bitget(v, 3));
    map_base{16} = @(v)(bitget(v, 2));
    map_base{17} = @(v)(bitget(v, 1));
elseif strcmp(name, 'F17xor')
    map_base = cell(25, 1);
    map_base{1} = @(v)(1);
    map_base{2} = @(v)(bitget(v, 16));
    map_base{3} = @(v)(bitget(v, 15));
    map_base{4} = @(v)(bitget(v, 14));
    map_base{5} = @(v)(bitget(v, 13));
    map_base{6} = @(v)(bitget(v, 12));
    map_base{7} = @(v)(bitget(v, 11));
    map_base{8} = @(v)(bitget(v, 10));
    map_base{9} = @(v)(bitget(v, 9));
    map_base{10} = @(v)(bitget(v, 8));
    map_base{11} = @(v)(bitget(v, 7));
    map_base{12} = @(v)(bitget(v, 6));
    map_base{13} = @(v)(bitget(v, 5));
    map_base{14} = @(v)(bitget(v, 4));
    map_base{15} = @(v)(bitget(v, 3));
    map_base{16} = @(v)(bitget(v, 2));
    map_base{17} = @(v)(bitget(v, 1));
    
    map_base{18} = @(v)(bitxor(map_base{2}(v), map_base{10}(v)));
    map_base{19} = @(v)(bitxor(map_base{3}(v), map_base{11}(v)));
    map_base{20} = @(v)(bitxor(map_base{4}(v), map_base{12}(v)));
    map_base{21} = @(v)(bitxor(map_base{5}(v), map_base{13}(v)));
    map_base{22} = @(v)(bitxor(map_base{6}(v), map_base{14}(v)));
    map_base{23} = @(v)(bitxor(map_base{7}(v), map_base{15}(v)));
    map_base{24} = @(v)(bitxor(map_base{8}(v), map_base{16}(v)));
    map_base{25} = @(v)(bitxor(map_base{9}(v), map_base{17}(v)));
elseif strcmp(name, 'F17tran')
    map_base = cell(33, 1);
    map_base{1} = @(v)(1);
    map_base{2} = @(v)(bitget(v, 16));
    map_base{3} = @(v)(bitget(v, 15));
    map_base{4} = @(v)(bitget(v, 14));
    map_base{5} = @(v)(bitget(v, 13));
    map_base{6} = @(v)(bitget(v, 12));
    map_base{7} = @(v)(bitget(v, 11));
    map_base{8} = @(v)(bitget(v, 10));
    map_base{9} = @(v)(bitget(v, 9));
    map_base{10} = @(v)(bitget(v, 8));
    map_base{11} = @(v)(bitget(v, 7));
    map_base{12} = @(v)(bitget(v, 6));
    map_base{13} = @(v)(bitget(v, 5));
    map_base{14} = @(v)(bitget(v, 4));
    map_base{15} = @(v)(bitget(v, 3));
    map_base{16} = @(v)(bitget(v, 2));
    map_base{17} = @(v)(bitget(v, 1));
    
    map_base{18} = @(v)(bitand(~map_base{2}(v), logical(map_base{10}(v))));
    map_base{19} = @(v)(bitand(~map_base{3}(v), logical(map_base{11}(v))));
    map_base{20} = @(v)(bitand(~map_base{4}(v), logical(map_base{12}(v))));
    map_base{21} = @(v)(bitand(~map_base{5}(v), logical(map_base{13}(v))));
    map_base{22} = @(v)(bitand(~map_base{6}(v), logical(map_base{14}(v))));
    map_base{23} = @(v)(bitand(~map_base{7}(v), logical(map_base{15}(v))));
    map_base{24} = @(v)(bitand(~map_base{8}(v), logical(map_base{16}(v))));
    map_base{25} = @(v)(bitand(~map_base{9}(v), logical(map_base{17}(v))));    
    
    map_base{26} = @(v)(bitand(logical(map_base{2}(v)), ~map_base{10}(v)));
    map_base{27} = @(v)(bitand(logical(map_base{3}(v)), ~map_base{11}(v)));
    map_base{28} = @(v)(bitand(logical(map_base{4}(v)), ~map_base{12}(v)));
    map_base{29} = @(v)(bitand(logical(map_base{5}(v)), ~map_base{13}(v)));
    map_base{30} = @(v)(bitand(logical(map_base{6}(v)), ~map_base{14}(v)));
    map_base{31} = @(v)(bitand(logical(map_base{7}(v)), ~map_base{15}(v)));
    map_base{32} = @(v)(bitand(logical(map_base{8}(v)), ~map_base{16}(v)));
    map_base{33} = @(v)(bitand(logical(map_base{9}(v)), ~map_base{17}(v)));  
else
    error('Unknown base name: %s', name);
end
    
