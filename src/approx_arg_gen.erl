%%%-------------------------------------------------------------------
%%% @author hrami
%%% @copyright (C) 2023, <COMPANY>
%%% @doc
%%%
%%% @end
%%% Created : 10. янв. 2023 19:10
%%%-------------------------------------------------------------------
-module(approx_arg_gen).

%% API
-export([generate_arguments/3]).

generate_arguments(Start, End, Step) ->
    generate_arguments([], Start, End, Step).

generate_arguments(Points, Current, End, Step) when Current =< End ->
    generate_arguments([Current | Points], Current + Step, End, Step);
generate_arguments(Points, _, _, _) ->
    Points.
