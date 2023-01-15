# Аппроксимация функций
Министерство науки и высшего образования Российской Федерации федеральное государственное автономное образовательное учреждение высшего образования

«Национальный исследовательский университет ИТМО»

---
__ФПИиКТ, Системное и Прикладное Программное Обеспечение__

__Лабораторная работа №3__

по Функциональному программированию

Выполнил: Хузин Р.Р.

Группа: P34112

Преподаватель: Пенской Александр Владимирович

###### Санкт-Петербург
###### 2022 г.
---

## Описание проблемы

Цель: получить навыки работы с вводом/выводом, потоковой обработкой данных, замыканиями, командной строкой.


Требования:

1. Программа должна быть реализована в функциональном стиле.
2. Ввод/вывод должен быть отделён от алгоритмов аппроксимации.
3. Алгоритм аппроксимации должен порождать замыкания, работающие на некоторых интервалах.
4. Требуется использовать идиоматичный для технологии стиль программирования.
5. Реализовать алгоритмы линейной, экспоненциальной и степенной аппроксимаций
6. Выбор пользователем необходимых алгоритмов
7. Настройка пользователем генератора вычисляемых точек


## Ключевые элементы реализации с минимальными комментариями

1. __Генератор вычисляемых точек__
    ```erlang
    generate_arguments(Start, End, Step) ->
      generate_arguments([], Start, End, Step).
    
    generate_arguments(Points, Current, End, Step) when Current =< End ->
        generate_arguments([Current | Points], Current + Step, End, Step);
    generate_arguments(Points, _, _, _) ->
        Points.
    ```
2. __Линейная аппроксимация__
   ```erlang
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
   ```
3. __Экспоненциальная аппроксимация__
    ```erlang
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
    ```

4. __Степенная аппроксимация__
    ```erlang
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
    ```

## Ввод программы
* __Код__
    ```erlang
    process_points(Points, Params) ->
    case io:fread("", "~f,~f") of
        {ok, [X, Y]} ->
            NewPoints = lists:append(Points, [#point{x = X, y = Y}]),
            case length(NewPoints) > 1 of
                true ->
                    [OldLastPoint | _] = lists:reverse(Points),
                    case OldLastPoint#point.x < X of
                        true ->
                            XValues =
                                approx_arg_gen:generate_arguments(OldLastPoint#point.x,
                                                                  X,
                                                                  Params#approx_params.step),
                            run_approx_funcs(Params#approx_params.approx_funcs,
                                             NewPoints,
                                             XValues,
                                             Params),
                            io:fwrite("----------------------------------~n"),
                            process_points(NewPoints, Params);
                        _ ->
                            io:fwrite("A next points X value must be greater than previous~n"),
                            process_points(NewPoints, Params)
                    end;
                _ ->
                    process_points(NewPoints, Params)
            end;
        eof ->
            ok;
        {error, Reason} ->
            exit(Reason)
    end.
    ```
* __Пример ввода данных__
    ```
     type src/input.txt | _build/default/bin/approx_funcs -step 0.3 -funcs lep -mp 4
  ```
  
    Файл input.txt:
    ```
    2.1,3.3
    3.4,1.2
    4.5,90.0
    5.321,41.2
    6.77,88.2
    7.1,89.1
  ```

## Вывод программы
* __Код__
    ```erlang
    write_results(Points, ApproxStr) ->
    io:fwrite("~s: ~n", [ApproxStr]),
    io:fwrite("----------------------~n"),
    write_points(Points),
    io:fwrite("----------------------~n~n~n").

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
    ```
* __Пример вывода данных__
    ```
    Power approximation: 
    ----------------------
    2.100000,3.300000
    2.400000,2.493228
    2.700000,1.947015
    3.000000,1.560641
    3.300000,1.277617
    ----------------------
    
    Exponential approximation:
    ----------------------
    2.100000,3.300000
    2.400000,2.612940
    2.700000,2.068926
    3.000000,1.638176
    3.300000,1.297108
    ----------------------
    
    Linear approximation:
    ----------------------
    2.100000,3.300000
    2.400000,2.815385
    2.700000,2.330769
    3.000000,1.846154
    3.300000,1.361538
    ----------------------
    
    ----------------------------------
    Power approximation:
    ----------------------
    3.400000,9.051973
    3.700000,12.310829
    4.000000,16.346107
    4.300000,21.263393
    ----------------------
    
    Exponential approximation:
    ----------------------
    3.400000,7.738441
    3.700000,11.473279
    4.000000,17.010677
    4.300000,25.220615
    ----------------------
    
    Linear approximation:
    ----------------------
    3.400000,33.832794
    3.700000,44.330370
    4.000000,54.827945
    4.300000,65.325520
    ----------------------
    
    ----------------------------------
    Power approximation:
    ----------------------
    4.500000,23.919511
    4.800000,30.075128
    5.100000,37.293339
    ----------------------
    
    Exponential approximation:
    ----------------------
    4.500000,22.938844
    4.800000,31.871259
    5.100000,44.281967
    ----------------------
    
    Linear approximation:
    ----------------------
    4.500000,47.145362
    4.800000,53.067136
    5.100000,58.988911
    ----------------------
    
    ----------------------------------
    Power approximation:
    ----------------------
    5.321000,42.930700
    5.621000,58.883336
    5.921000,79.447466
    6.221000,105.617699
    6.521000,138.537821
    ----------------------
    
    Exponential approximation:
    ----------------------
    5.321000,35.639225
    5.621000,49.473187
    5.921000,68.677034
    6.221000,95.335175
    6.521000,132.341120
    ----------------------
    
    Linear approximation:
    ----------------------
    5.321000,61.721778
    5.621000,67.820877
    5.921000,73.919975
    6.221000,80.019073
    6.521000,86.118172
    ----------------------
    
    ----------------------------------
    Power approximation:
    ----------------------
    6.770000,78.978703
    7.070000,80.641399
    ----------------------
    
    Exponential approximation:
    ----------------------
    6.770000,79.851071
    7.070000,82.240218
    ----------------------
    
    Linear approximation:
    ----------------------
    6.770000,82.178385
    7.070000,83.967721
    ----------------------
    
    ----------------------------------

  ```