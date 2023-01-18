-module(approx_funcs).

%% API exports
-export([main/1]).

-include("approx_funcs_methods.hrl").

-record(approx_params, {step = 0.2, approx_funcs = [l, e, p], points_window = 3}).

%%====================================================================
%% API functions
%%====================================================================

%% escript Entry point
main(Args) ->
    Params = parse_parameters(#approx_params{}, Args),
    Subscribers = spawn_generators(Params),
    process_points(Params, Subscribers),
    erlang:halt(0).

parse_parameters(Result, [Arg, NextArg | Tail]) ->
    case Arg of
        "-step" ->
            case string:to_float(NextArg) of
                {error, _} ->
                    io:fwrite("Wrong -step argument value: ~p ~n", [NextArg]),
                    exit(1);
                {Float, _} ->
                    parse_parameters(Result#approx_params{step = Float}, Tail)
            end;
        "-pw" ->
            case string:to_float(NextArg) of
                {error, _} ->
                    io:fwrite("Wrong -pw argument value: ~p ~n", [NextArg]),
                    exit(1);
                {Integer, _} ->
                    parse_parameters(Result#approx_params{points_window = Integer}, Tail)
            end;
        "-funcs" ->
            WithLinearFunc =
                case string:find(NextArg, "l") of
                    nomatch ->
                        Result#approx_params{approx_funcs = []};
                    _ ->
                        Result#approx_params{approx_funcs = [l]}
                end,
            WithExpFunc =
                case string:find(NextArg, "e") of
                    nomatch ->
                        WithLinearFunc;
                    _ ->
                        WithLinearFunc#approx_params{approx_funcs =
                                                         [e
                                                          | WithLinearFunc#approx_params.approx_funcs]}
                end,
            WithPowerFunc =
                case string:find(NextArg, "p") of
                    nomatch ->
                        WithExpFunc;
                    _ ->
                        WithExpFunc#approx_params{approx_funcs =
                                                      [p | WithExpFunc#approx_params.approx_funcs]}
                end,
            parse_parameters(WithPowerFunc, Tail)
    end;
parse_parameters(Result, [_ | _]) ->
    Result;
parse_parameters(Result, []) ->
    Result.


process_points(Params, Subscribers) ->
    process_points({0, true}, Params, Subscribers).

process_points({PreviousXValue, FirstRead}, Params, Subscribers) ->
    case io:fread("", "~f,~f") of
        {ok, [X, Y]} ->
            case (PreviousXValue < X) or FirstRead of
                true ->
                    send_point(#point{x = X, y = Y}, Subscribers),
                    wait_and_process_replies(Subscribers),
                    process_points({X, false}, Params, Subscribers);
                _ ->
                    io:fwrite("A next points X value must be greater than previous~n"),
                    process_points({PreviousXValue, false}, Params, Subscribers)
            end;
        eof ->
            ok;
        {error, Reason} ->
            exit(Reason)
    end.

send_point(Point, Subscribers) ->
    lists:foreach(
        fun({_, Sub}) ->
            Sub ! {send_point, self(), Point}
        end,
        Subscribers
    ).

wait_and_process_replies(Subscribers) ->
    lists:foreach(
        fun({Tag, _}) ->
            receive
                window_not_filled -> ok;
                not_possible -> write_approx_not_possible_result(Tag);
                {result, Points} -> write_approx_results(Points, Tag)
            after 5000 -> exit(timeout)
            end
        end,
        Subscribers
    ).

spawn_generators(Params) ->
    lists:map(
        fun(Atom) ->
            case Atom of
                l -> {l, spawn_link(
                    fun() -> approx_generator(fun approx_funcs_methods:calc_linear_least_squares/1, [], Params) end)
                };
                e -> {e, spawn_link(
                    fun() -> approx_generator(fun approx_funcs_methods:calc_exp_least_squares/1, [], Params) end)
                };
                p -> {p, spawn_link(
                    fun() -> approx_generator(fun approx_funcs_methods:calc_power_least_squares/1, [], Params) end)
                }
            end
        end,
        Params#approx_params.approx_funcs
    ).

approx_generator(Generator, PointsWindow, Params) ->
    receive
        {send_point, SenderPid, Point} ->
            NewWindow =
                case length(PointsWindow) > Params#approx_params.points_window - 1 of
                    true ->
                        [_ | Tail] = PointsWindow,
                        lists:append(Tail, [Point]);
                    _ ->
                        lists:append(PointsWindow, [Point])
            end,
            case length(NewWindow) >= Params#approx_params.points_window of
                true ->
                    [OldLastPoint | _] = lists:reverse(PointsWindow),
                    XValues =
                        approx_arg_gen:generate_arguments(OldLastPoint#point.x,
                            Point#point.x,
                            Params#approx_params.step),
                    Interval = {lists:nth(1, XValues), lists:last(XValues)},
                    case approx_funcs_methods:build_on_interval(NewWindow, Generator, Interval, Params#approx_params.points_window) of
                        not_possible ->
                            SenderPid ! not_possible;
                        Func ->
                            Result = approx_funcs_methods:calc_arguments(Func, XValues),
                            SenderPid ! {result, Result}
                    end;
                _ ->
                    SenderPid ! window_not_filled
            end,
            approx_generator(Generator, NewWindow, Params)
    end.


write_approx_results(Points, ApproxTag) ->
    case ApproxTag of
        l -> io:fwrite("~s: ~n", ["Linear approximation"]);
        e -> io:fwrite("~s: ~n", ["Exponential approximation"]);
        p -> io:fwrite("~s: ~n", ["Power approximation"])
    end,
    io:fwrite("----------------------~n"),
    write_points(Points),
    io:fwrite("----------------------~n~n").

write_approx_not_possible_result(ApproxTag) ->
    case ApproxTag of
        l -> io:fwrite("Linear approximation is not possible~n~n");
        e -> io:fwrite("Exponential approximation is not possible~n~n");
        p -> io:fwrite("Power approximation is not possible~n~n")
    end.

write_points([]) ->
    ok;
write_points([Point | Tail]) ->
    #point{x = X, y = Y} = Point,
    case Y of
        undef ->
            io:fwrite("~f,undef~n", [X]);
        _ ->
            io:fwrite("~f,~f~n", [X, Y])
    end,
    write_points(Tail).

%%====================================================================
%% Internal functions
%%====================================================================
