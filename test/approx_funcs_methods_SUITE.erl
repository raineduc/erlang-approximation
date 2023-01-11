%%%-------------------------------------------------------------------
%%% @author hrami
%%% @copyright (C) 2023, <COMPANY>
%%% @doc
%%%
%%% @end
%%% Created : 11. янв. 2023 9:58
%%%-------------------------------------------------------------------
-module(approx_funcs_methods_SUITE).

-include("approx_funcs_methods.hrl").

-include_lib("common_test/include/ct.hrl").
-include_lib("stdlib/include/assert.hrl").

%% API
-export([all/0]).
-export([test_linear_approx/1, test_exp_approx/1, test_power_approx/1, test_not_possible_cases/1]).

all() -> [test_linear_approx, test_exp_approx, test_power_approx, test_not_possible_cases].

%%--------------------------------------------------------------------
%% TEST CASES
%%--------------------------------------------------------------------

test_linear_approx(_) ->
  Points = [
    #point{x = 3.1, y = -1.9},
    #point{x = -2.2, y = 7.76},
    #point{x = 51.6, y = 10.151}
  ],
  LinearFunc = approx_funcs_methods:calc_linear_least_squares(Points),
  ?assert(abs(LinearFunc(-2.5) - 2.82725) < 1.0e-3),
  ?assert(abs(LinearFunc(-1) - 3.0155) < 1.0e-3),
  ?assert(abs(LinearFunc(4.43) - 3.696965) < 1.0e-3).

test_exp_approx(_) ->
  Points = [
    #point{x = 3.1, y = 1.9},
    #point{x = -2.2, y = 7.76},
    #point{x = 51.6, y = 10.151}
  ],
  ExpFunc = approx_funcs_methods:calc_exp_least_squares(Points),
  ?assert(abs(ExpFunc(-2.5) - 3.7994) < 1.0e-3),
  ?assert(abs(ExpFunc(-1) - 3.8958) < 1.0e-3),
  ?assert(abs(ExpFunc(4.43) - 4.26559) < 1.0e-3).

test_power_approx(_) ->
  Points = [
    #point{x = 3.1, y = 1.9},
    #point{x = 2.2, y = 7.76},
    #point{x = 51.6, y = 10.151}
  ],
  PowerFunc = approx_funcs_methods:calc_power_least_squares(Points),
  ?assertEqual(PowerFunc(-2.5), undef),
  ?assertEqual(PowerFunc(-0.1), undef),
  ?assert(abs(PowerFunc(1.1) - 3.14081) < 1.0e-3),
  ?assert(abs(PowerFunc(3.56) - 4.37604) < 1.0e-3),
  ?assert(abs(PowerFunc(137.512) - 12.2806) < 1.0e-2).

test_not_possible_cases(_) ->
  Points = [
    #point{x = -3.1, y = 1.9},
    #point{x = 2.2, y = -7.76},
    #point{x = 51.6, y = 10.151}
  ],
  PowerFunc = approx_funcs_methods:calc_power_least_squares(Points),
  ExpFunc = approx_funcs_methods:calc_exp_least_squares(Points),
  ?assertEqual(PowerFunc, not_possible),
  ?assertEqual(ExpFunc, not_possible).
