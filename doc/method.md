# Функции

## Прандтля-Майера

$$
\sqrt{\frac{\gamma+1}{\gamma-1}} \cdot \arctan{\left( \frac{(\gamma-1) \cdot \left(M_{\mu}^2-1\right)}{\gamma+1}  \right)} - \arctan{\left( M_{\mu}^2-1 \right)}
$$


$$
q(\lambda) = {\left( {\gamma_{\lambda}+1} \over 2 \right)}^{\dfrac{1}{\gamma_{\lambda}-1}} \cdot \dfrac{\sqrt{\dfrac{\gamma_{\lambda}+1}{2}} \cdot M_{\lambda}}{\sqrt{1+ \dfrac{\gamma_{\lambda}-1}{2}\cdot M_{\lambda}^2}}\ \cdot {\left( 1 - \dfrac{\gamma_{\lambda}-1}{\gamma_{\lambda}+1} \cdot {\left( \sqrt{\dfrac{\gamma_{\lambda}+1}{2}} \cdot M_{\lambda} \over \sqrt{1+ \dfrac{\gamma_{\lambda}-1}{2}\cdot M_{\lambda}^2} \right)}^2 \right)}^{\dfrac{1}{\gamma_{\lambda}-1}} 
$$

## Число Маха в критике

$$
\alpha = \sqrt{\dfrac{2}{(\gamma+1) \cdot \rho}},
$$

$$
\nu_s = \varepsilon_s = \dfrac{1}{8} \sqrt{2 (\gamma+1)\dfrac{1}{\rho}},
$$

$$
M_{min} = \sqrt{1 + {(\alpha \cdot \nu_s)}^2},
$$

где $\rho \approx 5$ - радиус переходной части между критикой и сз профилем.

$$
\alpha_0 = \arcsin{\dfrac{1}{M_{min}}}
$$

$$
cst = -\alpha_0 -
	\sqrt{\dfrac{\gamma+1}{\gamma-1}} \cdot
	\arctan{\left(
			\sqrt{\dfrac{\gamma-1}{\gamma+1}} \cdot
			\dfrac{1}{\tan{\alpha_0}} 
		\right)},
$$

$$
\theta(\alpha) = 
	\sqrt{\dfrac{\gamma+1}{\gamma-1}} \cdot
	\arctan{\left(
			\sqrt{\dfrac{\gamma-1}{\gamma+1}} \cdot
			\dfrac{1}{\tan{\alpha}}
		\right)} +
	\alpha +
	cst,
$$

## Вычисление Throat region Mach = 1+ curve 

| ![Красная кривая](img/region-mach=1.png) |
| :--: |
| *Рисунок 1: Красным выделен фронт потока с $M=1$* |

В приведённом коде происходит настройка геометрии изопериметрической поверхности в области сопла Лаваля, где число Маха $M = 1$ в горловине и плавно возрастает вдоль оси сопла.

### Разбиение радиуса на сегменты

Переменная `nChar` обозначает количество характеристик (линий). Массив `radii` хранит радиальные координаты точек на горловине, где $M = 1$.
 

Инициализация массива:

```python
radii = np.zeros(nChar)
radii[0] = 1  # Радиус в первой точке
radii[nChar - 1] = 0  # Радиус в последней точке
```

### Параметры сегментов кривой

Здесь задаются параметры кривой:

- `sn_r1` и `sn_r2` — параметры, которые разделяют горловину на две части.
- `qr` — коэффициент, контролирующий экспоненциальный характер уменьшения радиуса.
- `b1_r1` и `b1_r2` — коэффициенты для вычисления радиусов в каждой части кривой.

Формулы для коэффициентов:
$$
b1\_r1 = sn\_r1 \cdot \frac{qr - 1}{qr^{\left(\frac{nChar}{2} - 1\right)} - 1}
$$
$$
b1\_r2 = sn\_r2 \cdot \frac{qr - 1}{qr^{\left(\frac{nChar}{2}\right)} - 1}
$$

### Расчёт радиусов

1. **Уменьшение радиуса от 1 до `sn_r1`:**

    - Радиус уменьшается экспоненциально:
    ```python
    for i in range(1, int(nChar / 2)):
        r1_current = r1_current - b1_r1 * mt.pow(qr, i - 1)
        radii[i] = r1_current
    ```

	- Формула для радиуса:
    $$
    r_{1}(i) = r_{1\_current} - b1\_r1 \cdot qr^{(i-1)}
	$$

2. **Увеличение радиуса от `sn_r1` до 0:**

    - Радиус увеличивается экспоненциально:
    ```python
    for i in range(1, int(nChar / 2)):
        r2_current = r2_current + b1_r2 * mt.pow(qr, i - 1)
        radii[nChar - (i + 1)] = r2_current
    ```

    - Формула для радиуса:
    $$
    r_{2}(i) = r_{2\_current} + b1\_r2 \cdot qr^{(i-1)}
    $$

### Заполнение массива точек

В конце, массив `ArrayPtsIsoMach` заполняется следующими значениями:
- Радиальная координата `r`: \(\text{{radii[i]}}\)
- Координата вдоль оси `z`:
  $$z = 2 \cdot \epsilon_s \cdot \left(1 - r^2\right)$$
- Угол $\alpha$:
  $$\alpha = \arcsin{\left(\frac{1}{M\_min}\right)}$$
- Угол $\theta$ задается нулем.

```python
for i in range(nChar):
    ArrayPtsIsoMach[0, i, 0] = radii[i]  # r
    ArrayPtsIsoMach[0, i, 1] = 2 * eps_s * (1 - mt.pow(radii[i], 2))  # z
    ArrayPtsIsoMach[0, i, 2] = mt.asin(1 / M_min)  # alpha
    ArrayPtsIsoMach[0, i, 3] = 0  # theta
```

Таким образом, этот код формирует кривую для горловины сопла, где число Маха в каждой точке равно или больше 1, и описывает геометрию соответствующего изопериметрического сечения.