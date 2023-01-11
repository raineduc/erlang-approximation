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
    read_points(Points) ->
      case io:fread("", "~f,~f") of
          {ok, [X, Y]} ->
              read_points([#point{x = X, y = Y} | Points]);
          eof ->
              Points;
          {error, Reason} ->
              {error, Reason}
      end.
    ```
* __Пример ввода данных__
    ```
    type src/input.txt | _build/default/bin/erlang_approx_funcs -is -3.0 -ie -1.0 -step -funcs e
  ```
  
    Файл input.txt:
    ```
    2.1,3.3
    4.4,1.2
    -1.1,90.0
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
    Power approximation is not possible
    
    Exponential approximation:
    ----------------------
    -3.000000,330.824987
    -2.900000,305.371472
    -2.800000,281.876338
    -2.700000,260.188908
    -2.600000,240.170099
    -2.500000,221.691527
    -2.400000,204.634688
    -2.300000,188.890193
    -2.200000,174.357073
    -2.100000,160.942124
    -2.000000,148.559314
    -1.900000,137.129232
    -1.800000,126.578575
    -1.700000,116.839681
    -1.600000,107.850092
    -1.500000,99.552158
    -1.400000,91.892663
    -1.300000,84.822485
    -1.200000,78.296284
    -1.100000,72.272206
    ----------------------
    
    
    Linear approximation:
    -2.300000,100.520249
    -2.200000,98.836828
    -2.100000,97.153408
    -2.000000,95.469987
    -1.900000,93.786566
    -1.800000,92.103145
    -1.700000,90.419725
    -1.600000,88.736304
    -1.500000,87.052883
    -1.400000,85.369463
    -1.300000,83.686042
    -1.200000,82.002621
    -1.100000,80.319201
    ----------------------
  ```