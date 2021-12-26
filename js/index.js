// скрипт программы расчета струйного компрессора v1.1
function q(lum,k){ //функция приведенной массовой скорости;
    var rMV = Math.pow(((k+1)/2),(1/(k-1)))*lum*Math.pow((1-(k-1)/(k+1)*lum*lum),(1/(k-1)));
    return rMV;
}
function rP(lum,k){ //функция приведенного давления;
    var reducedPressure = Math.pow((1-(k-1)/(k+1)*lum*lum),(k/(k-1)));
    return reducedPressure;
}
function rIV(P,k){ //функция приведенной изоэнтропной скорости;
    var reducedIsentropicVelocity = Math.sqrt((1+k-Math.pow(P,(k-1)/k)*(k+1))/(k-1));
    return reducedIsentropicVelocity;
}
function lum1(q1,k){ //функция приведенной изоэнтропной скорости в диапазоне дозвуковых скоростей;
    for(zi=0,q=0;zi<=100000;zi++){
        var li = zi/100000;
        var q = Math.pow((k+1)/2,1/(k-1))*li*Math.pow((1-(k-1)/(k+1)*Math.pow(li,2)),1/(k-1));
        if (Math.abs((q-q1)/q1)<0.0001){return li;
             break;
            //  continue
            };
        if (zi==100000){return 1}
         }
}
function lum2(q1,k){ //функция приведенной изоэнтропной скорости в диапазоне сверхзвуковых скоростей;
    //var promtr = Math.sqrt((k+1)/(k-1))*100000;
    for(yi=100000,q=0;yi<Math.sqrt((k+1)/(k-1))*100000;yi++){
        var li = yi/100000;
        var q = Math.pow((k+1)/2,1/(k-1))*li*Math.pow(1-((k-1)/(k+1)*li*li),1/(k-1));
         if (Math.abs((q-q1)/q1)<0.0001){
                                            return li;
                                         }
         }
}
function diam(square){ //функция нахождения диаметра (в мм) по известной площади поперечного сечения канала (в м.кв);
    var d = Math.sqrt(4*square/Math.PI)*Math.pow(10, 3);
    return d;
}
function ro(p,t,r){ //функция расчета плотности газа по уравнению Касперовича (давление в Па, температура в К);
    var delta = molecularWeightGas.value /(22.4 * 1.293);
    var z = 1-5.5*(p*Math.pow(delta,1.3)/Math.pow(t,3.3));
    var roG = p/(z*r*t);
    return roG;
    }
function vkp(kk,p,r){//функция определения критической скорости
   var t = 2*(kk/(kk+1))*(p/r);
   var kp = Math.sqrt(t);
   return kp;
}
function getMaxOfArray(numArray) {
    return Math.max.apply(null, numArray);
  }

function coloring(x){ //функция раскраски влияющих элементов на применение методики (она таки не заработала и ею надо пожертвовать) 
    if(x==0){                                    
        results.style.backgroundColor = "red";
        injectionPressure.style.backgroundColor = "red";
        mixturePressure.style.backgroundColor = "red";}
    else {
        results.style.backgroundColor = "";
        injectionPressure.style.backgroundColor = "";
        mixturePressure.style.backgroundColor = "";}
    return results,injectionPressure,mixturePressure;
}
function un(y1,y2,y3,k1,k2,r1,r2,a1kr,a2kr,pPs,qPN,φ3,lumPN,kp1,kp2,t1,t2, step){//функция основного расчета
    var arrLam=[]; var arrQc3=[];  var arrP=[]; var arrU=[]; var arrQn2=[];
    var arrLamH2=[]; var arrK3=[]; var arrAckr=[]; var arrU2=[]; var arrmin=[];
    var arrR3=[]; var arrT3=[]; var arrro3=[];
    var qc3=0; var upr; var u=1; var qn2=0; var kp3=0; var kp4=0; var u2=1; 
    var jj;
    if (step==0.1){jj=10};
    if (step==0.05){jj=20};
    if (step==0.01){jj=100};
    for(var j=1;j<=jj;j++){
        la = j/jj
        upr=0.4;
        var z=0;
        arrLam[j-1]=la;
        while(Math.abs((upr-Math.min(u,u2))/Math.min(u,u2))>0.01){
            z=z+1;
            upr = Math.min(u,u2);
            var k3 = ((k1/(k1-1))+(upr*(k2/(k2-1))*(r2/r1)))/((1/(k1-1))+(upr*(1/(k2-1))*(r2/r1)));
            qc3 = q(la,k3);
            var ackr = Math.sqrt((Math.pow(a1kr,2)+u*Math.pow(a2kr,2))/(1+u));
            var lynps = rIV(pPs,k3);
            var qps = q(lynps,k3);
            var qc3ned = y1/y3*qps;
            if(z>55) break;
            if(qc3>qc3ned) break;
            if(qc3>=y2/y3){
                u = (ackr/a2kr*k2/k3*rP(1,k2)/rP(1,k3)*y2/y3/qc3-a1kr/a2kr*k2/k1*rP(1,k2)/rP(1,k1)*y2/y1/qps)/(1-(ackr/a2kr*k2/k3*rP(1,k2)/rP(1,k3)*y2/y3/qc3));
                qn2 = u/((1+u)*(ackr/a2kr*k2/k3*rP(1,k2)/rP(1,k3)*y2/y3/qc3)-(a1kr/a2kr*k2/k1*rP(1,k2)/rP(1,k1)*y2/y1*1/qPN));
                kp3=1+φ3*y3/y1*a1kr/ackr*(rP(lum1(qc3,k3),k3)-y2/y3)/(k1*rP(1,k1)*lum1(qc3,k3)*q(lumPN,k1));
                kp4=1+φ3*a2kr/ackr*y3/y2*(rP(lum1(qc3,k3),k3)-rP(lum1(qn2,k2),k2)*y2/y3)/(k2*rP(1,k2)*lum1(qc3,k3)*qn2);
                u2=(kp1*a1kr/ackr*lumPN-kp3*lum1(qc3,k3))/(kp4*lum1(qc3,k3)-kp2*a2kr/ackr*lum1(qn2,k3));
                if(u>u2){
                    u=0;
                    z1=0;
                    while (Math.abs((u-u2)/u2)>0.01){
                        u=u2;
                        z1=z1+1;
                        qn2 = u/((1+u)*(ackr/a2kr*k2/k3*rP(1,k2)/rP(1,k3)*y2/y3/qc3)-(a1kr/a2kr*k2/k1*rP(1,k2)/rP(1,k1)*y2/y1*1/qPN));
                        kp3=1+φ3*y3/y1*a1kr/ackr*(rP(lum1(qc3,k3),k3)-y2/y3)/(k1*rP(1,k1)*lum1(qc3,k3)*q(lumPN,k1));
                        kp4=1+φ3*a2kr/ackr*y3/y2*(rP(lum1(qc3,k3),k3)-rP(lum1(qn2,k2),k2)*y2/y3)/(k2*rP(1,k2)*lum1(qc3,k3)*qn2);
                        u2=(kp1*a1kr/ackr*lumPN-kp3*lum1(qc3,k3))/(kp4*lum1(qc3,k3)-kp2*a2kr/ackr*lum1(qn2,k3));
                        if ((u<0.001) && (u2<0.001)) break;
                        if(z1>55) break;
                        }
                    }
                 }
             }
        arrQc3[j-1]=+qc3.toFixed(3);
        arrP[j-1]=+rP(qc3,k3).toFixed(3);
        arrU[j-1]=+u.toFixed(3);
        arrQn2[j-1]=+qn2.toFixed(3);
        arrLamH2[j-1]=+lum1(qn2,k2).toFixed(3);
        arrK3[j-1]=+k3.toFixed(3);
        arrAckr[j-1]=+ackr.toFixed(3);
        arrU2[j-1]=+u2.toFixed(3);
        arrmin[j-1]=+Math.min(u,u2).toFixed(3);
        var r3 = (r1+arrmin[j-1]*r2)/(1+arrmin[j-1]);
        arrR3[j-1]=+r3.toFixed(3);
        var t3 = (t1+arrmin[j-1]*t2)/(1+arrmin[j-1]);
        arrT3[j-1]=+t3.toFixed(3);
        var ro3=+ro(y3,t3,r3);
        arrro3[j-1]=+ro3.toFixed(3)
    }
    var resultIt = {
        1: arrLam,
        2: arrQc3,
        3: arrP,
        4: arrU,
        5: arrQn2,
        6: arrLamH2,
        7: arrK3,
        8: arrAckr,
        9: arrU2,
        10: arrmin,
        11: arrR3,
        12: arrT3,
        13: arrro3
        };
        var resultKays = Object.keys(resultIt);
        for (y2=0;y2<resultIt[resultKays[3]].length;y2++){
            if ((resultIt[resultKays[3]][y2]==1)||(resultIt[resultKays[3]][y2]<=0.001)){
                resultIt[resultKays[3]][y2]=0;
                resultIt[resultKays[4]][y2]=0;
                resultIt[resultKays[5]][y2]=0;
                resultIt[resultKays[6]][y2]=0;
                resultIt[resultKays[7]][y2]=0;
                resultIt[resultKays[8]][y2]=0;
                resultIt[resultKays[9]][y2]=0;
                resultIt[resultKays[10]][y2]=r1;
                resultIt[resultKays[11]][y2]=t1;
                resultIt[resultKays[12]][y2]=+ro(y3,t1,r1).toFixed(3);
                                        }
           }
    return resultIt;
    }
function vuvod(res,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15){// функция определения столбца содержащего umax;
                //x1=g1,x2=a1kr,x3=k1,x4=p1,x5=ppn,x6=ro1,x7=v1,x8=p3,x9=qPN,x10=k2,x11=a2kr,
                //x12=p2,x13=β, x14=v3, x15=r
    var maxi=Math.max(...res[10]);
    for (rr=0;rr<res[10].length;rr++){
        if (res[10][rr]==maxi)break;
        }
    var λc3= res[1][rr];
    var qc3 = res[2][rr];
    var qн2 = res[5][rr];
    var k3 = res[7][rr];
    var aскр = res[8][rr];
    var u = res[10][rr];
    var ro3= res[13][rr];

    var fpKp = x1*x2/(x3*rP(1,x3)*x4);//площадь критического сечения сопла, м.кв;
    var qp1 = q(rIV(x5,x3),x3); //приведенная массовая скорость эжектирующего потока в выходном сечении сопла;
    var fp1 = fpKp/qp1;// площадь выходного сечения  сопла, м.кв;
    var fp = x1/(x6*x7);// площадь входного сечения сопла, м.кв;
    var f3 = (x3/k3*rP(1,x3)/rP(1,k3)*x4*aскр*(1+u)/(x8*x2*qc3))*fpKp; //площадь поперечного сечения камеры смешения
    var f2 = fpKp*(1/x9+x3/x10*rP(1,x3)/rP(1,x10)*(x4*x11*u)/(x12*x2*qн2));//площадь входного сечения камеры смешения, м.кв(для проверки f2=f3);
    var g1 = fp1/fpKp;//коэффициент подобия сопла;
    var g2 = f3/fpKp;//коэффициент подобия сопла;
    var dkr = diam(fpKp);//диаметр критического сечения сопла, мм;
    var dp1 = diam(fp1);//диаметр выходного сечения сопла, мм;
    var dp = diam(fp);//диаметр входного сечения сопла, мм;
    var d3 = diam(f3);//диаметр выходного сечения камеры смешения, мм;
    ///расчет геометрических осевых размеров;
    var a//опытная константа;
            if (u<0.2){a=0.07}
            else{
                if (0.2<=u<0.5){a=0.08}
                else {a=0.09};
            };
    var lc1//длинна свободной струи, мм;
    var d4//диаметр свободной струи на расстоянии lc1, мм;
            if (u<0.5){
                        let ff=Math.pow((0.083+0.76*u),0.5);
                        lc1=(ff-0.29)*dp1/(2*a);
                        d4=3.4*dp1*ff
                        }
            else {
                    lc1=(0.37+u)/(4.4*a)*dp1;
                    d4=1.55*dp1*(1+u);
                };
    var dkam = Math.max(d3,d4)//диаметр камеры смешения, мм
    var lc //Расстояние выходного сечения сопла от входного сечения цилинтдрической камеры смешения, мм
            if(d3>d4){lc=lc1}
            else {lc=lc1+(d4-d3)/(2*Math.tan(Math.PI*x13/180))};
    var lk = 8*d3; //длина цилиндрической камеры смешения для струйных аппаратов 
                    //(выбирается в пределах 6-10 диаметров камеры смешения), мм
    var fc = (x1*(1+u))/(ro3*x14);//площадь на выходе из диффузора, м.кв
    var dc = diam(fc);//диаметр выходного сечения диффузора, мм
    var ld = 7*(dc-d3);//длина диффузора, исходя из угла раствора 8-10 градусов (коэффициент = 6-7), мм
    var q2 = x1*u;//массовый расход инжектируемого потока, кг/сек
    var v2 = q2/x15*24*60*60; // объемный расход инжектируемого потока, кг/сек

    var maxiMal = {
        n: ++rr, λc3: λc3, qc3: qc3, qн2: qн2, k3: k3, aс_кр: aскр, u: u, fp_кр: fpKp, qp1: qp1.toFixed(3),
        fp1: fp1, fp: fp, f3:f3, f2:f2, g1: g1.toFixed(3), g2:g2.toFixed(3),
        dkr: dkr.toFixed(1), dp1: dp1.toFixed(1), dp: dp, d3: d3.toFixed(1), a: a, lc1: lc1.toFixed(1), d4: d4.toFixed(1),
        dкам: dkam.toFixed(1), lc: lc.toFixed(1), lk: lk.toFixed(1), fc: fc.toFixed(1),
        dc: dc.toFixed(1), ld: ld.toFixed(1), Q2: q2.toFixed(0), V2: v2.toFixed(0)
        }
    return maxiMal
}
var styleWidth = ["5%", "65%", "10%", "10%", "10%"]; //массив значений ширин столбцов таблиц
//функция отрисовки таблицы результатов расчета струйного компрессора;
function resu(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27,x28,
                x29,x30,x31,x32,x33){
    const resultsDate = document.getElementById('resultsDate')
    resultsDate.classList.remove('resultsDate')
    var parameters = {
        0: ["Параметр", "Единицы измерения", "Обозначение", "Значение"],
        1: ["Степень сжатия эжектируеющего потока","-","К<sub>с</sub>",x1.toFixed(3)],
        2: ["Степень расширения инжектируемого потока","-","К<sub>p</sub>",x2.toFixed(3)],
        3: ["Плотность эжектируеющего потока","кг/м<sup>3</sup>", "r<sub>1</sub>",x3.toFixed(3)],
        4: ["Плотность инжектируемого потока","кг/м<sup>3</sup>", "r<sub>2</sub>",x4.toFixed(3)],
        5: ["Объемный расход эжектируеющего потока при рабочих условиях","м<sup>3</sup>/с", "Q1",x5.toFixed(3)],
        6: ["Массовый расход эжектируеющего потока при рабочих условиях","м<sup>3</sup>/с", "G1",x6.toFixed(3)],
        7: ["Скорость рабочего потока в подводящем трубопроводе","м/с", "V<sub>1</sub>",x7.toFixed(3)],
        8: ["Критическая скорость эжектируеющего потока","м/с", "a<sub>1кр</sub>",x8.toFixed(3)],
        9: ["Критическая скорость инжектируемого потока","м/с", "a<sub>2кр</sub>",x9.toFixed(3)],
        10: ["Относительная критическая скорость","-", "Θ",x10.toFixed(3)],
        11: ["Относительное давление","-", "Π<sub>.р.н.</sub>",x11.toFixed(3)],
        12: ["Приведенная изоэнтропная скорость эжектируеющего потока","-", "λ<sub>.р.н.</sub>",x12.toFixed(3)],
        13: ["Приведенная массовая скорость эжектируеющего потока","-", "q<sub>.р.н.</sub>",x13.toFixed(3)],
        14: ["Приведенное давление в критическом сечении","-", "П<sub>н.кр</sub>",x14.toFixed(3)],
        15: ["Приведенное давление в любом сечении цил. камеры смешения","-", "Π<sub>р.s.</sub>",x15.toFixed(3)],
        16: ["Коэффициент скорости эжектируеющего потока","-", "K<sub>.1</sub>",x16.toFixed(3)],
        17: ["Коэффициент скорости инжектируемого потока","-", "K<sub>.2</sub>",x17.toFixed(3)],
        18: ["Порядковый номер столбца, в котором находится максимальное значение коэффициента эжекции","-", "n",x18],
        19: ["Оптимальная приведенная массовая скорость сжатого потока в сечении 3-3","-", "q<sub>.с3</sub>",x19],
        20: ["Оптимальная приведенная массовая скорость сжатого потока в сечении 2-2","-", "q<sub>.н2</sub>",x20],
        21: ["Максимальный коэффициент эжекции","-", "U",x21],
        22: ["Объемный расход инжектируемого потока","н.м<sup>3</sup>/сут", "Q<sub>2.сут</sub>",x22],
        23: ["Критический диаметр сопла","мм","d<sub>.kr</sub>",x23],
        24: ["Диаметр выходного сечения сопла","мм","d<sub>p1</sub>",x24],
        25: ["Диаметр камеры смешения","мм","d<sub>kam</sub>",x25],
        26: ["Диаметр выходного сечения диффузора","мм","d<sub>c</sub>",x26],
        27: ["Расстояние выходного сечения сопла от входного сечения цилинтдрической камеры смешения","мм","l<sub>.c</sub>",x27],
        28: ["Длина цилиндрической камеры","мм","l<sub>.k</sub>",x28],
        29: ["Длина диффузора","мм","l<sub>.д</sub>",x29],
        30: ["Время начала расчета","чч:мм:сс:мс","t<sub>нач</sub>",x30],
        32: ["Время затраченное на расчет","мс","t<sub>зат</sub>",x32],
        31: ["Время завершения расчета","чч:мм:сс:мс","t<sub>зав</sub>",x31],
        33: ["Версия программы расчета","","","v.1.01"],
        34: ["Детальный отчет", "", "", "<button id = 'detalies' class = 'calc' >Просмотреть</button>"]
                        }
    var forecast = Object.keys(parameters);
    //var results = document.createElement("div");
    var calculationResults = document.createElement("table");
    var zag = document.createElement("caption");
        zag.innerHTML=" 5. Результаты расчета";
        calculationResults.appendChild(zag);
    for (y=0;y<(forecast.length);y++){
        if (y>0){
            var tr = document.createElement("tr");
            var td = document.createElement("td");
            td.innerHTML = forecast[y]+".";
//            td.onmousemove = function(){
//            td.parentElement.style.backgroundColor = "blue";
//                }
            tr.appendChild(td);
            for (i=0;i<parameters[y].length;i++){
                                                    if (parameters[y][i] !== undefined) {
                                                        var td = document.createElement("td");
                                                        td.innerHTML = parameters[y][i];
                                                        //td.onmousemove = function(){ //подсвет строки ;-)
                                                        //    td.parentElement.style.backgroundColor = "green";
                                                        //    };
                                                        //td.onmouseout = function(){ //отбой подсветки
                                                        //    td.style.backgroundColor = ""
                                                        //    };
                                                        if (i==0){td.style.textAlign= "left"};
                                                        tr.appendChild(td)
                                                        }
                                                    else {break};
                                                };
                    }
        else{
            var tr = document.createElement("tr");
            var th = document.createElement("th");
            th.innerHTML = "№ п./п.";
            th.style.width = styleWidth[y];
            tr.appendChild(th);
            for (i=0;i<parameters[y].length;i++){
                    var th = document.createElement("th");
                    th.innerHTML = parameters[y][i];
                    th.style.width = styleWidth[i+1];
                    tr.appendChild(th)
                    };
            };
            calculationResults.appendChild(tr);
                                    };
    results.appendChild(calculationResults);
    //return results;
detalies.onclick = function(){//функция вывода развернутого отчета по кнопке 
    var detaliesResults = document.createElement("table");
    //detaliesResults.innerHTML = "" //чистим таблицу с расчетом;
    var zag = document.createElement("caption");
        zag.innerHTML=" 6. Результаты расчета развернуто";
        detaliesResults.appendChild(zag);
        x33["1"].unshift('λ<sub>c3</sub>');
        x33["2"].unshift('q<sub>c3</sub>(<sub>λc3</sub>)');
        x33["3"].unshift('П<sub>c3</sub>');
        x33["4"].unshift('U<sub>пр2</sub>');
        x33["5"].unshift('q<sub>н2</sub>');
        x33["6"].unshift('λ<sub>н2</sub>');
        x33["7"].unshift('k<sub>3</sub>');
        x33["8"].unshift('a<sub>с_кр</sub>');
        x33["9"].unshift('U<sub>2</sub>');
        x33["10"].unshift('U');
        x33["11"].unshift('R<sub>3</sub>');
        x33["12"].unshift('T<sub>3</sub>');
        x33["13"].unshift('r<sub>3</sub>');
        var a = x33["1"].length;
        var b = Object.keys(x33).length;
        for (y=0;y<b;y++){
            if (y==0){
                var tr = document.createElement("tr");
                for (i=0;i<a;i++){
                    var th = document.createElement("td");
                    th.innerHTML = x33["1"][i];
                    tr.appendChild(th);
                    };
                detaliesResults.appendChild(tr);
                }
            else{
                var tr = document.createElement("tr");
                for (i=0;i<a;i++){
                    var td = document.createElement("td");
                    td.innerHTML = x33[y+1][i];
                    tr.appendChild(td);
                    };
                detaliesResults.appendChild(tr);
                };
        }
        results.appendChild(detaliesResults);
        var p = document.createElement("p");
        results.appendChild(p);
    };
    return calculationResults;  
}

calc.onclick = function(){
    var time1 = new Date ()// переменная фиксирующая время начала расчета;
    var time11 = time1.getHours()+" :" + time1.getMinutes() +" :" + time1.getSeconds() +" :" + time1.getMilliseconds()
    results.innerHTML = "" //чистим таблицу с расчетом;
    var m0 = +molecularWeightGas.value; // Получаем M0;
    var p1 = +workflofPressure.value*Math.pow(10, 6);  //получаем ввод Р1 и переводим из МПа в Па;
    var t1 = +workflofTemperature.value; //получаем T1;
    var q1n = +massFlowWorkflow.value; //массовый расход эжектирующего потока;
    var k1 = +indicatorAdiabatEjectingFlow.value; //показатель адиабаты эжектирующего потока;
    var r1 = +unevenGasConstantEjectingFlow.value;//получаем R1;
    var r = +gasDensityUnderNC.value; //плотность рабочего потока при нормальных условиях;
    var p2 = +injectionPressure.value*Math.pow(10, 6);  //получаем ввод Р2 и переводим из МПа в Па;
    var t2 = +injectionFlowTemperature.value; //получаем T2;
    var k2 = +indicatorAdiabatInjectedStream.value; //показатель адиабаты инжектируюемого потока;
    var r2 = +unevenGasConstantInjectedStream.value;//получаем R2;
    var p3 = +mixturePressure.value*Math.pow(10, 6); //получаем ввод Р3 и переводим из МПа в Па;
    var v3 = +mixSpeed.value// Скорость смеси на выходе из камеры смешения, м/с;
    var step = stepλc3.value// шаг итераций λc3;
    var dp = +dpod.value; // получаем диаметр трубопровода подводящего эжектирующий поток;
    var β = +beta.value;//Угол между образующей входного участка камеры смешения и осью компрессора;
    var φ1 = +fi1.value; //Коэффициент скорости сопла
    var φ2 = +fi2.value; //Коэффициент скорости камеры смещения
    var φ3 = +fi3.value; //Коэффициент скорости диффузора
    var φ4 = +fi4.value; //	Коэффициент скорости конфузора
    //Предварительные рассчеты
    var kc = p3/p2;//Степень сжатия инжектируемого потока
    var kp = p1/p3;//Степень расширения эжектируемого потока
    var ro1 = ro(p1,t1,r1);//плотность эжектирующего потока 
    var ro2 = ro(p2,t2,r2);//плотность инжектируемого потока 
    var q1 = q1n*r/ro1/24/60/60;//Объемный расход эжектирующего потока  при рабочих условиях
    var g1 = q1n*r/ro1/24/60/60*ro1;//Массовый расход инжектируемого потока при рабочих условиях
    var v1 = q1/(Math.PI*(Math.pow(dp*Math.pow(10,-3),2))/4);//скорость эжектирующего потока, м/с
    var a1kr = vkp(k1,p1,ro1);//критическая скорость эжектирующего потока
    var a2kr = vkp(k2,p2,ro2);//критическая скорость инжектируемого потока
    var tet = Math.pow((a2kr/a1kr),2);//Относительная критическая скорость
    var ppn = p2/p1;//Относительное давление
    var lumPN = rIV(ppn,k1);//приведенная изоэнтропная скорость эжектирующего потока
    var qPN = q(lumPN,k1);//приведенная массовая скорость эжектирующего потока
    var pNkr = rP(1,k2);//приведенное давление в критическом сечении
    var pPs = ppn*pNkr;// приведенное давление в любом сечении цил. камеры смешения
    var kp1 = φ1*φ2*φ3;//Коэффициент скорости эжектирующего потока
    var kp2 = φ1*φ2*φ4;//Коэффициент скорости инжектируемого потока

    // проверка применимости методики расчёта основной расчет и вывод результатов в таблицу
    var p11 = document.createElement("p");
    if ((kc.toFixed(1)<=2.5)&&(kc.toFixed(1)>=1.2)){
        p11.innerHTML = "Применение методики разрешено";
        var p12 = document.createElement("p");
        p12.innerHTML = "1,2 <= K<sub>c</sub> <= 2,5";
        results.appendChild(p11);
        results.appendChild(p12);
        var uu2=un(p1,p2,p3,k1,k2,r1,r2,a1kr,a2kr,pPs,qPN,φ3,lumPN,kp1,kp2,t1,t2,step);//вызов функции основного расчета
        //console.log(p1,p2,p3,uu2);
        var uI=vuvod(uu2,g1,a1kr,k1,p1,ppn,ro1,v1,p3,qPN,k2,a2kr,p2,β,v3,r);//вызов функции нахождения максимального значения
        //console.log(uI);
        var time2 = new Date ()////вставить переменную фиксирующую время начала расчета;
        var time21 = time2.getHours()+" :" + time2.getMinutes() +" :" + time2.getSeconds() +" :" + time2.getMilliseconds();
        var times = time2 - time1;
        //var times1 = times.getHours()+" :" + times.getMinutes() +" :" + times.getSeconds() +" :" + times.getMilliseconds();
        resu(kc,kp,ro1,ro2,q1,g1,v1,a1kr,a2kr,tet,ppn,lumPN,qPN,pNkr,pPs,kp1,kp2,uI.n,uI.qc3,uI.qн2,uI.u,uI.V2,uI.dkr,uI.dp1,uI.dкам,uI.dc,uI.lc,uI.lk,uI.ld,time11,time21,times,uu2);// отрисовываем таблицу результатов расчета
        var p22 = document.createElement("p");
        results.appendChild(p22);
        coloring(1); //функция подсветки влияющих ячеек и полей
       }
    else {if (kc.toFixed(1)>=2.5) {
                                    p11.style.backgroundColor = "red";
                                    p11.innerHTML = "Применима методика расчета эжектора."
                                    results.appendChild(p11);
                                    coloring(0); //функция  подсветки влияющих ячеек и полей
                                                }
            else {
                    p11.style.backgroundColor = "red";
                    p11.innerHTML = "Применима методика расчета инжектора."
                    results.appendChild(p11);
                    coloring(0);//функция подсветки влияющих ячеек и полей
                }
    };
}