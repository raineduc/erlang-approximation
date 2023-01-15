%%%-------------------------------------------------------------------
%%% @author hrami
%%% @copyright (C) 2023, <COMPANY>
%%% @doc
%%%
%%% @end
%%% Created : 10. янв. 2023 13:28
%%%-------------------------------------------------------------------
-module(approx_funcs_methods).

-include("approx_funcs_methods.hrl").

%% Partial result of linear approximation
%% sx - Sum of X values
%% sxx - sum of X to the power of 2 values
%% sy - sum of Y values
%% sxy - sum of X*Y values
-record(linap_partial, {sx, sxx, sy, sxy, number_of_points}).
%% Partial result of exponential approximation
%% sx - Sum of X values
%% sxx - sum of X to the power of 2 values
%% sy - sum of log(Y) values
%% sxy - sum of X*log(Y) values
-record(expap_partial, {sx, sxx, sy, sxy, number_of_points}).
%% Partial result of power approximation
%% sx - Sum of log(X) values
%% sxx - sum of log(X) to the power of 2 values
%% sy - sum of log(Y) values
%% sxy - sum of log(X)*log(Y) values
-record(powap_partial, {sx, sxx, sy, sxy, number_of_points}).

%% API
-export([calc_linear_least_squares/1, calc_exp_least_squares/1,
         calc_power_least_squares/1, calc_arguments/2, build_on_interval/4]).

calc_linear_least_squares(Points) when length(Points) > 1 ->
    PartialResult =
        #linap_partial{sx = 0,
                       sxx = 0,
                       sy = 0,
                       sxy = 0,
                       number_of_points = length(Points)},
    calc_linear_least_squares(Points, PartialResult).

calc_linear_least_squares([],
                          #linap_partial{sx = SX,
                                         sxx = SXX,
                                         sy = SY,
                                         sxy = SXY,
                                         number_of_points = NumberOfPoints}) ->
    Determinant = SXX * NumberOfPoints - SX * SX,
    case Determinant of
        0 ->
            not_possible;
        _ ->
            CoefA = (SXY * NumberOfPoints - SX * SY) / Determinant,
            CoefB = (SXX * SY - SX * SXY) / Determinant,
            fun(X) -> CoefA * X + CoefB end
    end;
calc_linear_least_squares([Point | Tail],
                          #linap_partial{sx = SX,
                                         sxx = SXX,
                                         sy = SY,
                                         sxy = SXY,
                                         number_of_points = NumberOfPoints}) ->
    #point{x = X, y = Y} = Point,
    calc_linear_least_squares(Tail,
                              #linap_partial{sx = SX + X,
                                             sxx = SXX + X * X,
                                             sy = SY + Y,
                                             sxy = SXY + X * Y,
                                             number_of_points = NumberOfPoints}).

calc_exp_least_squares(Points) when length(Points) > 1 ->
    PartialResult =
        #expap_partial{sx = 0,
                       sxx = 0,
                       sy = 0,
                       sxy = 0,
                       number_of_points = length(Points)},
    calc_exp_least_squares(Points, PartialResult).

calc_exp_least_squares([],
                       #expap_partial{sx = SX,
                                      sxx = SXX,
                                      sy = SY,
                                      sxy = SXY,
                                      number_of_points = NumberOfPoints}) ->
    Determinant = SXX * NumberOfPoints - SX * SX,
    case Determinant of
        0 ->
            not_possible;
        _ ->
            CoefA = math:exp((SXX * SY - SX * SXY) / Determinant),
            CoefB = (SXY * NumberOfPoints - SX * SY) / Determinant,
            fun(X) -> CoefA * math:exp(CoefB * X) end
    end;
calc_exp_least_squares([Point | Tail],
                       #expap_partial{sx = SX,
                                      sxx = SXX,
                                      sy = SY,
                                      sxy = SXY,
                                      number_of_points = NumberOfPoints}) ->
    #point{x = X, y = Y} = Point,
    case Y =< 0 of
        true ->
            not_possible;
        _ ->
            calc_exp_least_squares(Tail,
                                   #expap_partial{sx = SX + X,
                                                  sxx = SXX + X * X,
                                                  sy = SY + math:log(Y),
                                                  sxy = SXY + X * math:log(Y),
                                                  number_of_points = NumberOfPoints})
    end.

calc_power_least_squares(Points) when length(Points) > 1 ->
    PartialResult =
        #powap_partial{sx = 0,
                       sxx = 0,
                       sy = 0,
                       sxy = 0,
                       number_of_points = length(Points)},
    calc_power_least_squares(Points, PartialResult).

calc_power_least_squares([],
                         #powap_partial{sx = SX,
                                        sxx = SXX,
                                        sy = SY,
                                        sxy = SXY,
                                        number_of_points = NumberOfPoints}) ->
    Determinant = SXX * NumberOfPoints - SX * SX,
    case Determinant of
        0 ->
            not_possible;
        _ ->
            CoefA = math:exp((SXX * SY - SX * SXY) / Determinant),
            CoefB = (SXY * NumberOfPoints - SX * SY) / Determinant,
            fun(X) ->
               case X < 0 of
                   true ->
                       undef;
                   _ ->
                       CoefA * math:pow(X, CoefB)
               end
            end
    end;
calc_power_least_squares([Point | Tail],
                         #powap_partial{sx = SX,
                                        sxx = SXX,
                                        sy = SY,
                                        sxy = SXY,
                                        number_of_points = NumberOfPoints}) ->
    #point{x = X, y = Y} = Point,
    case {X =< 0, Y =< 0} of
        {true, _} ->
            not_possible;
        {_, true} ->
            not_possible;
        _ ->
            calc_power_least_squares(Tail,
                                     #powap_partial{sx = SX + math:log(X),
                                                    sxx = SXX + math:log(X) * math:log(X),
                                                    sy = SY + math:log(Y),
                                                    sxy = SXY + math:log(X) * math:log(Y),
                                                    number_of_points = NumberOfPoints})
    end.

build_on_interval(Points, ApproxBuilder, _, MaxPoints) when length(Points) =< MaxPoints ->
    ApproxBuilder(Points);
build_on_interval(Points, ApproxBuilder, {IntervalStart, IntervalEnd}, MaxPoints) ->
    ApproxWindow = choose_approx_window({1, MaxPoints}, Points, {IntervalStart, IntervalEnd}),
    ApproxBuilder(ApproxWindow).

choose_approx_window({WindowsStartIndex, WindowEndIndex}, Points, Interval)
    when WindowEndIndex < length(Points) ->
    CurrentMetric =
        calc_window_convergence_metric({WindowsStartIndex, WindowEndIndex}, Points, Interval),
    NextMetric =
        calc_window_convergence_metric({WindowsStartIndex + 1, WindowEndIndex + 1},
                                       Points,
                                       Interval),
    case NextMetric < CurrentMetric of
        true ->
            choose_approx_window({WindowsStartIndex + 1, WindowEndIndex + 1}, Points, Interval);
        _ ->
            lists:sublist(Points, WindowsStartIndex, WindowEndIndex - WindowsStartIndex + 1)
    end;
choose_approx_window({WindowsStartIndex, WindowEndIndex}, Points, _) ->
    lists:sublist(Points, WindowsStartIndex, WindowEndIndex - WindowsStartIndex + 1).

calc_window_convergence_metric({WindowsStartIndex, WindowEndIndex}, Points, Interval) ->
    point_deviation(lists:nth(WindowsStartIndex, Points), Interval)
    + point_deviation(lists:nth(WindowEndIndex, Points), Interval).

point_deviation(#point{x = X}, {IntervalStart, IntervalEnd}) ->
    Left = X < IntervalStart,
    Right = X > IntervalEnd,
    case true of
      Left -> IntervalStart - X;
      Right -> X - IntervalEnd;
      _ -> 0
    end.

calc_arguments(ApproxFunc, Arguments) ->
    calc_arguments([], ApproxFunc, Arguments).

calc_arguments(ResultPoints, _, []) ->
    ResultPoints;
calc_arguments(ResultPoints, ApproxFunc, [X | Tail]) ->
    calc_arguments([#point{x = X, y = ApproxFunc(X)} | ResultPoints], ApproxFunc, Tail).
