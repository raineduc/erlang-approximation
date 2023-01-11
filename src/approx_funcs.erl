-module(approx_funcs).

%% API exports
-export([main/1]).

-include("approx_funcs_methods.hrl").

-record(approx_params, {interval_start = -2.0, interval_end = 2.0, step = 0.2, approx_funcs = [l, e, p]}).

%%====================================================================
%% API functions
%%====================================================================

%% escript Entry point
main(Args) ->
  io:format("Args: ~p~n", [Args]),
  Params = parse_parameters(#approx_params{}, Args),
  case read_points([]) of
    {error, Reason} -> exit(Reason);
    Points ->
      XValues = approx_arg_gen:generate_arguments(Params#approx_params.interval_start, Params#approx_params.interval_end, Params#approx_params.step),
      run_approx_funcs(Params#approx_params.approx_funcs, Points, XValues)
  end,
  erlang:halt(0).

parse_parameters(Result, [Arg, NextArg | Tail]) ->
  case Arg of
    "-is" ->
        case string:to_float(NextArg) of
          {error, _} ->
            io:fwrite("Wrong -is argument value: ~p ~n", [NextArg]),
            exit(1);
          {Float, _} -> parse_parameters(Result#approx_params{interval_start = Float}, Tail)
        end;
    "-ie" ->
        case string:to_float(NextArg) of
          {error, _} ->
            io:fwrite("Wrong -ie argument value: ~p ~n", [NextArg]),
            exit(1);
          {Float, _} -> parse_parameters(Result#approx_params{interval_end = Float}, Tail)
        end;
    "-step" ->
      case string:to_float(NextArg) of
        {error, _} ->
          io:fwrite("Wrong -step argument value: ~p ~n", [NextArg]),
          exit(1);
        {Float, _} -> parse_parameters(Result#approx_params{step = Float}, Tail)
      end;
    "-funcs" ->
      WithLinearFunc = case string:find(NextArg, "l") of
                         nomatch -> Result#approx_params{approx_funcs = []};
                         _ -> Result#approx_params{approx_funcs = [l]}
                       end,
      WithExpFunc = case string:find(NextArg, "e") of
                      nomatch -> WithLinearFunc;
                      _ -> WithLinearFunc#approx_params{approx_funcs = [e | WithLinearFunc#approx_params.approx_funcs]}
                    end,
      WithPowerFunc = case string:find(NextArg, "p") of
                        nomatch -> WithExpFunc;
                        _ -> WithExpFunc#approx_params{approx_funcs = [p | WithExpFunc#approx_params.approx_funcs]}
                      end,
      parse_parameters(WithPowerFunc, Tail)
  end;
parse_parameters(Result, [_ | _]) -> Result;
parse_parameters(Result, []) -> Result.

read_points(Points) ->
  case io:fread("", "~f,~f") of
    {ok, [X, Y]} ->
       read_points([#point{x = X, y = Y} | Points]);
    eof ->
      Points;
    {error, Reason} -> {error, Reason}
  end.


run_approx_funcs([], _, _) -> ok;
run_approx_funcs([FuncLiteral | Tail], ApproxPoints, XValues) ->
  case FuncLiteral of
    l ->
      case approx_funcs_methods:calc_linear_least_squares(ApproxPoints) of
        not_possible -> io:fwrite("Linear approximation is not possible~n~n");
        LinearFunc -> write_results(approx_funcs_methods:calc_arguments(LinearFunc, XValues), "Linear approximation")
      end;
    e ->
      case approx_funcs_methods:calc_exp_least_squares(ApproxPoints) of
        not_possible -> io:fwrite("Exponential approximation is not possible~n~n");
        ExpFunc -> write_results(approx_funcs_methods:calc_arguments(ExpFunc, XValues), "Exponential approximation")
      end;
    p ->
      case approx_funcs_methods:calc_power_least_squares(ApproxPoints) of
        not_possible -> io:fwrite("Power approximation is not possible~n~n");
        PowerFunc -> write_results(approx_funcs_methods:calc_arguments(PowerFunc, XValues), "Power approximation")
      end
  end,
  run_approx_funcs(Tail, ApproxPoints, XValues).

write_results(Points, ApproxStr) ->
  io:fwrite("~s: ~n", [ApproxStr]),
  io:fwrite("----------------------~n"),
  write_points(Points),
  io:fwrite("----------------------~n~n~n").

write_points([]) -> ok;
write_points([Point | Tail]) ->
  #point{x = X, y = Y} = Point,
  case Y of
    undef -> io:fwrite("~f,undef~n", [X]);
    _ -> io:fwrite("~f,~f~n", [X, Y])
  end,
  write_points(Tail).


%%====================================================================
%% Internal functions
%%====================================================================
